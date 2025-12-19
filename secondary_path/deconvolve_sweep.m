function out = deconvolve_sweep(recorded, sweepSig, fs, params)
% deconvolve_sweep (修复版 v1.4)
% 修复点：使用能量比例代替绝对值比例，添加侧瓣抑制检查，修复peakReliability判定
%
% 输出字段:
%   h                : 截断后的时域IR（最终用于ANC）
%   delayCorr        : 互相关估计延迟 (样本)
%   startIdxGlobal   : 截取开始在原始 recorded 中的索引
%   peakIdx          : h 内主峰索引 (1-based)
%   peakIdxZeroBased : 0-based 索引
%   peakReliability  : 布尔，峰可信
%   peakAbsFrac      : 峰值区域绝对值比例
%   peakEnergyFrac   : 峰值区域能量比例
%   sideLobeSuppression : 侧瓣抑制比 (dB)
%   preEnergyFrac    : 峰前能量占总能量比例
%   snrEst           : 主瓣窗口 vs 尾部 SNR (dB)
%   noiseStd, noiseMAD, thresholdUsed
%   triggerIdx       : 累计能量触发位置
%   pkLocalGlobal    : 原始 h_full 中的峰位置
%   warnEarly        : 是否出现录音早于播放的情况
%   paramsUsed       : 参数回传

if nargin < 4 || isempty(params), params = struct(); end

regEps        = getP(params,'regEps',1e-4);
extraTail     = getP(params,'extraTail',4096);      % 仅用于 rec_aligned 长度
preDelayKeep  = getP(params,'preDelayKeep',256);
tailTotal     = getP(params,'tailTotal',4096);      % 最终输出 IR 长度
peakThreshDB  = getP(params,'peakThreshDB',12);
maxSearch     = getP(params,'maxSearch',7000);      % 必须 >= 物理延迟！
noiseWin      = getP(params,'noiseWin',400);
envSmoothWin  = getP(params,'envSmoothWin',8);
cumEnergyFrac = getP(params,'cumEnergyFrac',0.05);
minPeakFrac   = getP(params,'minPeakFrac',0.02);
snrBodyRadius = getP(params,'snrBodyRadius',96);
fftCorrEnable = getP(params,'fftCorrEnable',true);
debugMode     = getP(params,'debugMode',false);     % 新增调试模式

rec = recorded(:);
exc = sweepSig(:);
Nrec = length(rec);
Nexc = length(exc);

%% 互相关延迟估计
if fftCorrEnable
    NfftCorr = 2^nextpow2(Nrec + Nexc);
    RECf = fft(rec, NfftCorr);
    EXCf = fft(exc, NfftCorr);
    C = ifft(RECf .* conj(EXCf));
    C = [C(end-Nexc+2:end); C(1:Nrec)];
    lags = (-Nexc+2:Nrec)';
else
    [C,lags] = xcorr(rec, exc);
end
[~,imax] = max(C);  % 不取绝对值
delayCorr = lags(imax);

startIdx = delayCorr + 1;
warnEarly = startIdx < 1;
if startIdx < 1, startIdx = 1; end

segEnd = min(startIdx + Nexc + extraTail - 1, Nrec);
rec_aligned = rec(startIdx:segEnd);

%% 频域反卷积 (Wiener 正则)
Nfft = 2^nextpow2(length(rec_aligned) + Nexc - 1);
REC = fft(rec_aligned, Nfft);
EXC = fft(exc, Nfft);
magEXC2 = abs(EXC).^2;
Hf = (REC .* conj(EXC)) ./ (magEXC2 + regEps);
h_full = real(ifft(Hf));

% 🔧 修复点：不再用 extraTail 截断 h_full！
maxLengthToKeep = max(maxSearch + 2048, 16384);  % 至少 16k，确保大延迟场景不丢峰
if length(h_full) > maxLengthToKeep
    h_full = h_full(1:maxLengthToKeep);
end

% 去局部均值（使用前段噪声估计）
nw_pre = min(noiseWin, floor(length(h_full)/10));
if nw_pre < 10, nw_pre = min(100, length(h_full)); end
noiseBaseMean = mean(h_full(1:nw_pre));
h_full = h_full - noiseBaseMean;

%% 噪声统计 + 阈值
nw = min(noiseWin, length(h_full)-10);
if nw < 10, nw = min(100, length(h_full)); end
noiseSeg = h_full(1:nw);
noiseStd = std(noiseSeg);
noiseMAD = median(abs(noiseSeg - median(noiseSeg))) / 0.6745;
noiseBase = max(noiseStd, noiseMAD);
th = noiseBase * 10^(peakThreshDB/20);

searchEnd = min(length(h_full), maxSearch);
if searchEnd < 10, searchEnd = min(length(h_full), 100); end

cand = h_full(1:searchEnd);
envRaw = abs(cand);
envSm = movmean(envRaw, envSmoothWin);
totalE = sum(envSm);
if totalE < 1e-18
    % 极弱信号，直接返回默认结构
    out = default_output(fs, params, tailTotal);
    return;
end

cumE = cumsum(envSm);
triggerIdx = find(cumE >= totalE * cumEnergyFrac, 1, 'first');
if isempty(triggerIdx), triggerIdx = 1; end

pkLocal = find(envSm(triggerIdx:searchEnd) > th, 1, 'first');
if isempty(pkLocal)
    [~,pkLocal] = max(envSm(triggerIdx:searchEnd));
    pkLocal = pkLocal + triggerIdx - 1;
else
    pkLocal = pkLocal + triggerIdx - 1;
end

%% ✅ 修复：改进的peakReliability判定
peakWinEnd = min(pkLocal + envSmoothWin*16, length(h_full));
peakWin = h_full(pkLocal:peakWinEnd);

% 1. 计算两种比例
peakAbsFrac = sum(abs(peakWin)) / (sum(abs(h_full)) + 1e-12);
peakEnergyFrac = sum(peakWin.^2) / (sum(h_full.^2) + 1e-12);

% 2. 侧瓣抑制比检查
side_lobe_radius = min(100, floor(length(h_full)/4));
side_start = max(1, pkLocal - side_lobe_radius);
side_end = min(length(h_full), pkLocal + side_lobe_radius);
exclude_radius = min(20, side_lobe_radius/2);
exclude_start = max(1, pkLocal - exclude_radius);
exclude_end = min(length(h_full), pkLocal + exclude_radius);

% 构建侧瓣区域（排除主峰区域）
side_region = [];
for i = side_start:side_end
    if i < exclude_start || i > exclude_end
        side_region = [side_region, i];
    end
end

if ~isempty(side_region)
    main_peak_val = max(abs(h_full(exclude_start:exclude_end)));
    max_side_lobe = max(abs(h_full(side_region)));
    if max_side_lobe > 0
        side_lobe_suppression_db = 20*log10(main_peak_val/(max_side_lobe + 1e-12));
    else
        side_lobe_suppression_db = Inf;
    end
else
    side_lobe_suppression_db = Inf;
end

% 3. 综合判定条件
abs_ok = (peakAbsFrac >= minPeakFrac);
energy_ok = (peakEnergyFrac >= minPeakFrac * 0.5);  % 能量阈值减半
side_lobe_ok = (side_lobe_suppression_db > 3);  % 至少3dB抑制

% 最终可靠性：满足(绝对值或能量条件) AND 有基本侧瓣抑制
peakReliability = (abs_ok || energy_ok) && side_lobe_ok;

% 调试信息
if debugMode
    fprintf('  [DEB-deconv] pk=%d, absFrac=%.4f, energyFrac=%.4f, sideSupp=%.1f dB -> reliable=%d\n', ...
        pkLocal, peakAbsFrac, peakEnergyFrac, side_lobe_suppression_db, peakReliability);
end

%% 截窗保留前导（生成最终输出 IR）
winStartLocal = max(pkLocal - preDelayKeep, 1);
winStopLocal  = min(winStartLocal + tailTotal - 1, length(h_full));
h_out = h_full(winStartLocal:winStopLocal);

peakIdxFinal = pkLocal - winStartLocal + 1;
if peakIdxFinal < 1, peakIdxFinal = 1; end

% SNR 自适应估计
bodyL = max(1, peakIdxFinal - snrBodyRadius);
bodyR = min(length(h_out), peakIdxFinal + snrBodyRadius);
bodySlice = h_out(bodyL:bodyR);
tailSlice = h_out(max(length(h_out)-4*snrBodyRadius,1):end);
snrEst = 20*log10((rms(bodySlice)+1e-12)/(rms(tailSlice)+1e-12));

preEnergy = sum(abs(h_out(1:peakIdxFinal-1)));
postEnergy = sum(abs(h_out(peakIdxFinal:end)));
preEnergyFrac = preEnergy / (preEnergy + postEnergy + 1e-12);

out = struct();
out.h = h_out;
out.delayCorr = delayCorr;
out.startIdxGlobal = startIdx;
out.peakIdx = peakIdxFinal;
out.peakIdxZeroBased = peakIdxFinal - 1;
out.peakReliability = peakReliability;
out.peakAbsFrac = peakAbsFrac;
out.peakEnergyFrac = peakEnergyFrac;
out.sideLobeSuppression = side_lobe_suppression_db;
out.preEnergyFrac = preEnergyFrac;
out.snrEst = snrEst;
out.noiseStd = noiseStd;
out.noiseMAD = noiseMAD;
out.thresholdUsed = th;
out.triggerIdx = triggerIdx;
out.pkLocalGlobal = pkLocal;
out.warnEarly = warnEarly;
out.sampleRate = fs;
out.paramsUsed = params;
end

function val = getP(p, name, defaultVal)
    if isstruct(p) && isfield(p,name) && ~isempty(p.(name))
        val = p.(name);
    else
        val = defaultVal;
    end
end

function out = default_output(fs, params, tailTotal)
    % 安全兜底：当信号极弱或异常时返回合理默认值
    out.h = zeros(tailTotal, 1);
    out.delayCorr = 0;
    out.startIdxGlobal = 1;
    out.peakIdx = 1;
    out.peakIdxZeroBased = 0;
    out.peakReliability = false;
    out.peakAbsFrac = 0;
    out.peakEnergyFrac = 0;
    out.sideLobeSuppression = 0;
    out.preEnergyFrac = 0;
    out.snrEst = -Inf;
    out.noiseStd = 0;
    out.noiseMAD = 0;
    out.thresholdUsed = 0;
    out.triggerIdx = 1;
    out.pkLocalGlobal = 1;
    out.warnEarly = false;
    out.sampleRate = fs;
    out.paramsUsed = params;
end