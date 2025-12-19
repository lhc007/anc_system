function irAvgOut = average_ir(irRepList, doAlign, opts)
% average_ir 改进版 v2.3
%
% 输入:
%   irRepList: [Lh x numMics x reps] —— 冲激响应重复测量数据
%   doAlign  : 是否对齐主峰（逻辑标量）
%   opts     : 可选参数结构体
%       .peakSearchLen      —— 峰搜索最大索引（默认：Lh，即全长度）
%       .alignRefStrategy   —— 对齐参考策略：'median' (default), 'min', 'mode'
%       .excludeLowSNR      —— 是否剔除低SNR重复（默认 false）
%       .snrThresholdDB     —— SNR阈值（dB，默认 6）
%       .useMics            —— 用于峰检测的麦克风索引（默认全部）
%
% 输出结构 irAvgOut:
%   .irAvg              —— 平均/中值IR [Lh x numMics]
%   .peaks              —— 每次重复的综合主峰位置 [reps x 1]
%   .refPeak            —— 用作对齐参考的峰位置（标量）
%   .validMask          —— 哪些重复被保留 [reps x 1 logical]
%   .alignApplied       —— 是否执行了对齐
%   .snrEstGlobal       —— 全局SNR估计（dB）
%   .drift              —— 有效重复间主峰漂移（样本）
%   .peaksAllMics       —— 每次重复、每麦的峰位置 [reps x numMics]
%   .repSNR             —— 每次重复的SNR [reps x 1]

if nargin < 2, doAlign = true; end
if nargin < 3, opts = struct(); end

[Lh, numMics, reps] = size(irRepList);

% === 参数解析 ===
peakSearchLen    = getO(opts, 'peakSearchLen', Lh);           % ← 关键：默认搜索全长
alignRefStrategy = getO(opts, 'alignRefStrategy', 'median');
excludeLowSNR    = getO(opts, 'excludeLowSNR', false);
snrThresholdDB   = getO(opts, 'snrThresholdDB', 6);
useMics          = getO(opts, 'useMics', []);

if isempty(useMics)
    useMics = 1:numMics;
else
    useMics = useMics(useMics >= 1 & useMics <= numMics);
    if isempty(useMics), useMics = 1:numMics; end
end

% === 特殊处理：仅一次重复 ===
if reps == 1
    irSingle = irRepList(:, :, 1);
    env = mean(abs(irSingle(:, useMics)), 2);
    if ~isempty(env)
        [~, pk] = max(env);
    else
        pk = 1;
    end
    
    % 初始化输出
    irAvgOut.irAvg = irSingle;
    irAvgOut.peaks = pk;
    irAvgOut.refPeak = pk;
    irAvgOut.validMask = true;
    irAvgOut.alignApplied = false;
    irAvgOut.snrEstGlobal = NaN;
    irAvgOut.drift = 0;
    irAvgOut.peaksAllMics = zeros(1, numMics);
    irAvgOut.repSNR = NaN;

    % 填充每麦克风的峰位置
    for m = 1:numMics
        if Lh > 0
            [~, pkm] = max(abs(irSingle(:, m)));
            irAvgOut.peaksAllMics(1, m) = pkm;
        else
            irAvgOut.peaksAllMics(1, m) = 1;
        end
    end
    return;
end

% === 多重复处理 ===
peaks = zeros(reps, 1);
repSNR = zeros(reps, 1);
peaksAllMics = zeros(reps, numMics);

for r = 1:reps
    % 综合包络（仅使用指定麦克风）
    env = mean(abs(irRepList(:, useMics, r)), 2);
    searchEnd = min(peakSearchLen, length(env));
    if searchEnd < 1, searchEnd = 1; end
    [~, pk] = max(env(1:searchEnd));
    peaks(r) = pk;

    % SNR 估算：自适应窗口防越界
    winBody = min(128, floor(Lh * 0.1));  % 主瓣窗长 ≤ 总长10%
    winTail = min(512, floor(Lh * 0.2));  % 尾部窗长 ≤ 总长20%
    bodyStart = max(1, pk - winBody);
    bodyEnd   = min(Lh, pk + winBody);
    tailStart = max(1, Lh - winTail + 1);
    tailEnd   = Lh;

    body = irRepList(bodyStart:bodyEnd, useMics, r);
    tail = irRepList(tailStart:tailEnd, useMics, r);
    repSNR(r) = 20 * log10((rms(body(:)) + 1e-12) / (rms(tail(:)) + 1e-12));

    % 每麦克风单独找峰（全长度）
    for m = 1:numMics
        if Lh > 0
            [~, pkm] = max(abs(irRepList(:, m, r)));
            peaksAllMics(r, m) = pkm;
        else
            peaksAllMics(r, m) = 1;
        end
    end
end

% === 有效性筛选 ===
validMask = true(reps, 1);
if excludeLowSNR
    validMask = repSNR >= snrThresholdDB;
    if ~any(validMask)
        warning('[average_ir] 所有重复 SNR 低于 %.1f dB，保留全部', snrThresholdDB);
        validMask(:) = true; % 保底不剔除
    end
end

peaksValid = peaks(validMask);
drift = double(max(peaksValid) - min(peaksValid));

% === 确定参考峰 ===
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
if refPeak < 1, refPeak = 1; end
if refPeak > Lh, refPeak = Lh; end

% === 对齐或直接堆叠 ===
irStack = irRepList(:, :, validMask);
repsValid = sum(validMask);

if doAlign && drift > 0
    aligned = zeros(Lh, numMics, repsValid);
    pv = peaks(validMask);
    for k = 1:repsValid
        shift = pv(k) - refPeak;
        h_r = irStack(:, :, k);
        if shift > 0
            % 向右移：前面补零，后面截断
            if shift < Lh
                aligned(:, :, k) = [zeros(shift, numMics); h_r(1:end-shift, :)];
            else
                aligned(:, :, k) = zeros(Lh, numMics);
            end
        elseif shift < 0
            % 向左移：前面截断，后面补零
            s = -shift;
            if s < Lh
                aligned(:, :, k) = [h_r(s+1:end, :); zeros(s, numMics)];
            else
                aligned(:, :, k) = zeros(Lh, numMics);
            end
        else
            aligned(:, :, k) = h_r;
        end
    end
    irWork = aligned;
    alignApplied = true;
else
    irWork = irStack;
    alignApplied = false;
end

% === 最终平均（使用中值抗异常值）===
irAvg = median(irWork, 3);

% === 全局 SNR 估算 ===
envGlob = mean(abs(irAvg(:, useMics)), 2);
[~, pkGlob] = max(envGlob);
bodyStartG = max(1, pkGlob - winBody);
bodyEndG   = min(Lh, pkGlob + winBody);
tailStartG = max(1, Lh - winTail + 1);
tailEndG   = Lh;
bodyGlob = irAvg(bodyStartG:bodyEndG, useMics);
tailGlob = irAvg(tailStartG:tailEndG, useMics);
snrEstGlobal = 20 * log10((rms(bodyGlob(:)) + 1e-12) / (rms(tailGlob(:)) + 1e-12));

% === 构造输出 ===
irAvgOut.irAvg = irAvg;
irAvgOut.peaks = peaks;
irAvgOut.refPeak = refPeak;
irAvgOut.validMask = validMask;
irAvgOut.alignApplied = alignApplied;
irAvgOut.snrEstGlobal = snrEstGlobal;
irAvgOut.drift = drift;
irAvgOut.peaksAllMics = peaksAllMics;
irAvgOut.repSNR = repSNR;

% === 提示信息 ===
if snrEstGlobal < 6
    fprintf('[average_ir] 警告：全局 SNR = %.2f dB 偏低\n', snrEstGlobal);
else
    fprintf('[average_ir] 全局 SNR = %.2f dB\n', snrEstGlobal);
end
end

% -------------------------------------------------------------------------
function val = getO(s, name, defaultVal)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        val = s.(name);
    else
        val = defaultVal;
    end
end