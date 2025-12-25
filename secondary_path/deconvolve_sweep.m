function out = deconvolve_sweep(recorded, sweepSig, cfg)
% 彻底修复的反卷积函数 - 解决互相关失败问题

fs = cfg.fs;
regEps        = cfg.deconvRegEps;
preDelayKeep  = cfg.deconvPreDelayKeep;
tailTotal     = cfg.irMaxLen;
noiseWin      = cfg.deconvNoiseWin;
debugMode     = cfg.deconvDebugMode;

rec = recorded(:);
exc = sweepSig(:);
Nrec = length(rec);
Nexc = length(exc);

if debugMode
    fprintf('  [DEB-start] Nrec=%d, Nexc=%d, fs=%d Hz\n', Nrec, Nexc, fs);
end

%% ✅ 彻底修复：基于能量检测的延迟估计（不再依赖互相关）
% 1. 计算能量包络
energy_rec = movmean(rec.^2, 500);
energy_exc = movmean(exc.^2, 500);

% 2. 找到激励信号的起始点
exc_threshold = 0.1 * max(energy_exc);
exc_start = find(energy_exc > exc_threshold, 1, 'first');
if isempty(exc_start)
    exc_start = 1;
end

% 3. 在录音信号中搜索相似的能量模式
search_range = 1:min(Nrec, 10*Nexc);  % 搜索前10倍sweep长度
best_corr = -Inf;
best_lag = 0;

for lag = search_range
    end_idx = min(lag + Nexc - 1, Nrec);
    if end_idx - lag + 1 < Nexc * 0.5
        continue;  % 窗口太小
    end
    
    % 计算能量相关性
    rec_seg = energy_rec(lag:end_idx);
    exc_seg = energy_exc(1:length(rec_seg));
    
    corr_val = corr(rec_seg, exc_seg);
    if corr_val > best_corr
        best_corr = corr_val;
        best_lag = lag;
    end
end

if best_corr > 0.3
    delayEst = best_lag - exc_start;
else
    % 回退方案：使用简单的阈值检测
    rec_threshold = 0.05 * max(energy_rec);
    rec_start = find(energy_rec > rec_threshold, 1, 'first');
    if isempty(rec_start)
        rec_start = 1;
    end
    delayEst = rec_start - exc_start;
end

% 确保延迟为正
delayEst = max(1, delayEst);

if debugMode
    fprintf('  [DEB-delay] 能量相关延迟估计: %d (corr=%.3f)\n', delayEst, best_corr);
end

%% 反卷积
startIdx = delayEst;
segEnd = min(startIdx + Nexc + tailTotal - 1, Nrec);
rec_aligned = rec(startIdx:segEnd);

Nfft = 2^nextpow2(length(rec_aligned) + Nexc - 1);
REC = fft(rec_aligned, Nfft);
EXC = fft(exc, Nfft);
magEXC2 = abs(EXC).^2;

% 使用标准正则化
Hf = (REC .* conj(EXC)) ./ (magEXC2 + regEps);
h_full = real(ifft(Hf));

% 限制长度
maxLength = min(10*tailTotal, length(h_full));
h_full = h_full(1:maxLength);

% 去均值
noiseStart = min(noiseWin, floor(length(h_full)/10));
if noiseStart < 10
    noiseStart = min(100, length(h_full));
end
noiseMean = mean(h_full(1:noiseStart));
h_full = h_full - noiseMean;

%% 峰值检测
h_abs = abs(h_full);
h_energy = h_abs.^2;

% 简单但可靠的峰值检测
[peakVal, peakIdx] = max(h_abs);

% 计算噪声水平（使用前10%作为噪声估计）
noiseLen = min(floor(length(h_full)*0.1), 1000);
noiseLevel = std(h_full(1:noiseLen));

% SNR估计
if noiseLevel > 0
    snrEst = 20*log10(peakVal/(noiseLevel + 1e-12));
else
    snrEst = 30;
end

% 峰值可靠性判定
% 1. 峰值是否足够高
peak_high_enough = (peakVal > noiseLevel * 10);
% 2. 能量集中度
total_energy = sum(h_energy);
window_radius = min(100, floor(length(h_full)*0.05));
peak_window = h_energy(max(1,peakIdx-window_radius):min(length(h_full),peakIdx+window_radius));
energy_in_peak = sum(peak_window);
energy_concentration = energy_in_peak / (total_energy + 1e-12);
energy_concentrated = (energy_concentration > 0.1);

peakReliability = peak_high_enough && energy_concentrated;

if debugMode
    fprintf('  [DEB-peak] idx=%d, val=%.4f, noise=%.4f, SNR=%.1f dB, energy_conc=%.3f, reliable=%d\n', ...
        peakIdx, peakVal, noiseLevel, snrEst, energy_concentration, peakReliability);
end

%% 输出窗口
winStart = max(1, peakIdx - preDelayKeep);
winEnd = min(winStart + tailTotal - 1, length(h_full));

h_out = h_full(winStart:winEnd);

% 补零或截断
if length(h_out) < tailTotal
    h_out = [h_out; zeros(tailTotal - length(h_out), 1)];
elseif length(h_out) > tailTotal
    h_out = h_out(1:tailTotal);
end

% 调整峰值位置
peakIdxFinal = peakIdx - winStart + 1;
if peakIdxFinal < 1, peakIdxFinal = 1; end
if peakIdxFinal > tailTotal, peakIdxFinal = tailTotal; end

out = struct();
out.h = h_out;
out.delayEst = delayEst;
out.peakIdx = peakIdxFinal;
out.peakReliability = peakReliability;
out.snrEst = snrEst;
out.energyConcentration = energy_concentration;
out.sampleRate = fs;
end