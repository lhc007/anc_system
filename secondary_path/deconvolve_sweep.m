function out = deconvolve_sweep(recorded, sweepSig, fs, params)
% deconvolve_sweep (改进版 v1.2)
% 稳健反卷积 + 峰检测 + 可信度指标输出
%
% 输出字段:
%   h                : 截断后的时域IR
%   delayCorr        : 互相关估计延迟 (样本)
%   startIdxGlobal   : 截取开始在原始 recorded 中的索引
%   peakIdx          : h 内主峰索引 (1-based)
%   peakIdxZeroBased : 0-based 索引
%   peakReliability  : 布尔，峰可信
%   peakEnergyFrac   : 峰附近能量占总体能量比例
%   preEnergyFrac    : 峰前能量占总能量比例
%   snrEst           : 主瓣窗口 vs 尾部 SNR (dB)
%   noiseStd, noiseMAD, thresholdUsed
%   triggerIdx       : 累计能量触发位置
%   pkLocalGlobal    : 原始 h_full 中的峰位置
%   warnEarly        : 是否出现录音早于播放的情况
%   paramsUsed       : 参数回传
%
if nargin < 4 || isempty(params), params = struct(); end
regEps        = getP(params,'regEps',1e-4);
extraTail     = getP(params,'extraTail',4096);
preDelayKeep  = getP(params,'preDelayKeep',256);
tailTotal     = getP(params,'tailTotal',4096);
peakThreshDB  = getP(params,'peakThreshDB',12);
maxSearch     = getP(params,'maxSearch',7000);
noiseWin      = getP(params,'noiseWin',400);
envSmoothWin  = getP(params,'envSmoothWin',8);
cumEnergyFrac = getP(params,'cumEnergyFrac',0.05);
minPeakFrac   = getP(params,'minPeakFrac',0.02);
snrBodyRadius = getP(params,'snrBodyRadius',96);
fftCorrEnable = getP(params,'fftCorrEnable',true);

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
h_full = h_full(1:extraTail);

% 去局部均值
noiseBaseMean = mean(h_full(1:min(noiseWin, floor(extraTail/10))));
h_full = h_full - noiseBaseMean;

%% 噪声统计 + 阈值
nw = min(noiseWin, length(h_full)-10);
noiseSeg = h_full(1:nw);
noiseStd = std(noiseSeg);
noiseMAD = median(abs(noiseSeg - median(noiseSeg))) / 0.6745;
noiseBase = max(noiseStd, noiseMAD);
th = noiseBase * 10^(peakThreshDB/20);

searchEnd = min(length(h_full), maxSearch);
cand = h_full(1:searchEnd);

envRaw = abs(cand);
envSm = movmean(envRaw, envSmoothWin);
totalE = sum(envSm);
cumE = cumsum(envSm);
triggerIdx = find(cumE >= totalE * cumEnergyFrac, 1, 'first');
if isempty(triggerIdx), triggerIdx = 1; end

pkLocal = find(envSm(triggerIdx:searchEnd) > th, 1, 'first');
if isempty(pkLocal)
    [~,pkLocal] = max(envSm(triggerIdx:searchEnd));
end
pkLocal = pkLocal + triggerIdx - 1;

% 峰能量比例
peakWin = h_full(pkLocal:min(pkLocal+envSmoothWin*16,length(h_full)));
peakEnergyFrac = sum(abs(peakWin)) / (sum(abs(h_full))+1e-12);

peakReliability = (peakEnergyFrac >= minPeakFrac);

%% 截窗保留前导
winStartLocal = max(pkLocal - preDelayKeep, 1);
winStopLocal  = min(winStartLocal + tailTotal - 1, length(h_full));
h_out = h_full(winStartLocal:winStopLocal);

peakIdxFinal = pkLocal - winStartLocal + 1;
if peakIdxFinal < 1, peakIdxFinal = 1; end

% SNR 自适应
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
out.peakEnergyFrac = peakEnergyFrac;
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