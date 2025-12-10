function irAvgOut = average_ir(irRepList, doAlign, opts)
% average_ir 改进版 v2.2
% 输入:
%   irRepList: [Lh x numMics x reps]
%   doAlign  : 是否对齐主峰
%   opts 可选结构:
%       .peakSearchLen    默认 3000
%       .alignRefStrategy 'median'|'min'|'mode'
%       .excludeLowSNR    默认 false
%       .snrThresholdDB   默认 6
%       .useMics          指定用于峰检测的麦克风索引 (默认全部)
%
% 返回 irAvgOut 结构:
%   irAvgOut.irAvg
%   irAvgOut.peaks            (每重复主峰)
%   irAvgOut.refPeak
%   irAvgOut.validMask
%   irAvgOut.alignApplied
%   irAvgOut.snrEstGlobal
%   irAvgOut.drift            (峰漂移)
%   irAvgOut.peaksAllMics     (可选：每重复每麦克峰位置)
%
% 若只有 1 次重复，直接返回结构；保持与测量脚本兼容。
%
if nargin < 2, doAlign = true; end
if nargin < 3, opts = struct(); end

peakSearchLen    = getO(opts,'peakSearchLen',3000);
alignRefStrategy = getO(opts,'alignRefStrategy','median');
excludeLowSNR    = getO(opts,'excludeLowSNR',false);
snrThresholdDB   = getO(opts,'snrThresholdDB',6);
useMics          = getO(opts,'useMics',[]);

[Lh,numMics,reps] = size(irRepList);
if isempty(useMics)
    useMics = 1:numMics;
else
    useMics = useMics(useMics>=1 & useMics<=numMics);
    if isempty(useMics), useMics = 1:numMics; end
end

if reps == 1
    irAvg = irRepList(:,:,1);
    irAvgOut = struct('irAvg',irAvg,'peaks',1,'refPeak',1,'validMask',true, ...
        'alignApplied',false,'snrEstGlobal',NaN,'drift',0,'peaksAllMics',1);
    return;
end

% 获取每重复的综合峰包络与 SNR
peaks = zeros(reps,1);
repSNR = zeros(reps,1);
peaksAllMics = zeros(reps,numMics);

for r=1:reps
    % 综合包络（只使用选定麦克）
    env = mean(abs(irRepList(:,useMics,r)),2);
    searchLen = min(peakSearchLen,length(env));
    [~,pk] = max(env(1:searchLen));
    peaks(r) = pk;
    % 简单 SNR
    body = irRepList(max(1,pk-128):min(Lh,pk+128),useMics,r);
    tail = irRepList(max(Lh-512,1):Lh,useMics,r);
    repSNR(r) = 20*log10( (rms(body(:))+1e-12)/(rms(tail(:))+1e-12) );
    % 每麦克峰
    for m=1:numMics
        [~,pkm] = max(abs(irRepList(:,m,r)));
        peaksAllMics(r,m)=pkm;
    end
end

validMask = true(reps,1);
if excludeLowSNR
    validMask = repSNR >= snrThresholdDB;
    if ~any(validMask)
        validMask(:)=true; % 保底不全剔除
    end
end

peaksValid = peaks(validMask);
drift = max(peaksValid) - min(peaksValid);

switch lower(alignRefStrategy)
    case 'median'
        refPeak = median(peaksValid);
    case 'min'
        refPeak = min(peaksValid);
    case 'mode'
        refPeak = mode(round(peaksValid));
    otherwise
        refPeak = median(peaksValid);
end
refPeak = round(refPeak);

irStack = irRepList(:,:,validMask);
repsValid = sum(validMask);

if doAlign
    aligned = zeros(Lh,numMics,repsValid);
    pv = peaks(validMask);
    for k=1:repsValid
        shift = pv(k) - refPeak;
        h_r = irStack(:,:,k);
        if shift > 0
            if shift < Lh
                aligned(:,:,k) = [zeros(shift,numMics); h_r(1:end-shift,:)];
            else
                aligned(:,:,k) = zeros(Lh,numMics);
            end
        elseif shift < 0
            s = -shift;
            if s < Lh
                aligned(:,:,k) = [h_r(s+1:end,:); zeros(s,numMics)];
            else
                aligned(:,:,k) = zeros(Lh,numMics);
            end
        else
            aligned(:,:,k) = h_r;
        end
    end
    irWork = aligned;
    alignApplied = true;
else
    irWork = irStack;
    alignApplied = false;
end

irAvg = median(irWork,3);

% 全局 SNR 基于平均包络峰
[~,pkGlob] = max(mean(abs(irAvg(:,useMics)),2));
bodyGlob = irAvg(max(1,pkGlob-128):min(Lh,pkGlob+128),useMics);
tailGlob = irAvg(max(Lh-512,1):Lh,useMics);
snrEstGlobal = 20*log10( (rms(bodyGlob(:))+1e-12)/(rms(tailGlob(:))+1e-12) );

irAvgOut = struct();
irAvgOut.irAvg = irAvg;
irAvgOut.peaks = peaks;
irAvgOut.refPeak = refPeak;
irAvgOut.validMask = validMask;
irAvgOut.alignApplied = alignApplied;
irAvgOut.snrEstGlobal = snrEstGlobal;
irAvgOut.drift = drift;
irAvgOut.peaksAllMics = peaksAllMics;
irAvgOut.repSNR = repSNR;

% 提示
if snrEstGlobal < 6
    fprintf('[average_ir] 警告：全局 SNR=%.2f dB 偏低\n', snrEstGlobal);
else
    fprintf('[average_ir] 全局 SNR=%.2f dB\n', snrEstGlobal);
end
end

function val = getO(s,name,defaultVal)
if isfield(s,name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end