function result = deconv_industrial(recorded, sweep, cfg)
% 反卷积函数，使用频域维纳滤波
% 专门为ANC次级路径测量优化
% 输入:
%   recorded: 记录信号 (列向量)
%   sweep: 激励信号 (列向量)
%   cfg: 配置结构体
% 输出:
%   result: 包含ir, snr, peakPos的结构体

% 参数提取
fs = cfg.fs;
regEps = cfg.deconvRegEps;
irMaxLen = cfg.irMaxLen;

% 确保列向量
rec = recorded(:);
exc = sweep(:);
Nrec = length(rec);
Nexc = length(exc);

% ============== 频域反卷积 ==============
% 选择FFT长度，确保足够覆盖且为2的幂次
Nfft = 2^nextpow2(Nrec + Nexc - 1);
REC = fft(rec, Nfft);
EXC = fft(exc, Nfft);
magEXC2 = abs(EXC).^2;

% ============== 改进的正则化策略 ==============
% 自适应正则化：考虑激励信号的能量分布
maxEXC2 = max(magEXC2);
if maxEXC2 == 0
    maxEXC2 = 1; % 避免除零
end

% 基础正则化项
baseReg = regEps * maxEXC2;

% 频率相关正则化：在激励信号能量低的频段加强正则化
freqWeights = 1 ./ (sqrt(magEXC2 / maxEXC2) + 0.1); % 能量越低，权重越大
freqWeights = min(freqWeights, 10); % 限制最大权重

% 最终正则化参数
regEpsAdaptive = baseReg * freqWeights;

% ============== 维纳反卷积 ==============
Hf = (REC .* conj(EXC)) ./ (magEXC2 + regEpsAdaptive);
h_full = real(ifft(Hf));

% ============== 时域处理 ==============
% 截断到合理长度
maxKeep = min(irMaxLen * 2, length(h_full)); % 保留足够长度用于分析
h_full = h_full(1:maxKeep);

% 去直流偏置（使用前100样本的均值）
if length(h_full) >= 100
    dc_offset = mean(h_full(1:100));
else
    dc_offset = mean(h_full);
end
h_full = h_full - dc_offset;

% ============== 峰值检测与质量评估 ==============
h_abs = abs(h_full);

% 确定搜索范围
searchStart = max(1, cfg.minPhysDelaySamples);
searchEnd = min(length(h_full), cfg.maxPhysDelaySamples);
if searchEnd <= searchStart
    searchEnd = min(length(h_full), searchStart + 2000);
end

searchRange = searchStart:searchEnd;
if isempty(searchRange)
    searchRange = 1:length(h_full);
end

% 在搜索范围内寻找主峰
[peakValLocal, peakIdxLocal] = max(h_abs(searchRange));
if isempty(peakIdxLocal)
    peakIdx = searchStart;
else
    peakIdx = searchRange(peakIdxLocal);
end

% 增强的峰值验证
% 使用局部中值滤波后的信号进行峰值细化
localWin = 21; % 局部窗口大小
halfWin = floor(localWin/2);
refinedPeakIdx = peakIdx;

if peakIdx > halfWin && peakIdx <= length(h_abs) - halfWin
    localSegment = h_abs(peakIdx-halfWin:peakIdx+halfWin);
    % 使用抛物线插值提高精度
    [~, localMaxIdx] = max(localSegment);
    if localMaxIdx > 1 && localMaxIdx < length(localSegment)
        y1 = localSegment(localMaxIdx-1);
        y2 = localSegment(localMaxIdx);
        y3 = localSegment(localMaxIdx+1);
        delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
        refinedLocalIdx = localMaxIdx + delta;
        refinedPeakIdx = peakIdx - halfWin - 1 + refinedLocalIdx;
    end
end

peakIdx = round(refinedPeakIdx);

% ============== SNR计算 ==============
% 信号能量（主峰附近窗口）
winRadius = min(100, round(0.001 * fs)); % 1ms窗口或100样本
sigStart = max(1, peakIdx - winRadius);
sigEnd   = min(length(h_full), peakIdx + winRadius);
if sigEnd > sigStart
    signalSegment = h_full(sigStart:sigEnd);
    signalEnergy = sum(signalSegment.^2);
else
    signalEnergy = 0;
end

% 噪声能量（使用前导静音段和尾部静音段）
% 前导噪声段
preNoiseStart = 1;
preNoiseEnd = min(100, sigStart - 1);
if preNoiseEnd > preNoiseStart
    preNoise = h_full(preNoiseStart:preNoiseEnd);
else
    preNoise = [];
end

% 尾部噪声段（IR之后的区域）
postNoiseStart = min(sigEnd + 100, length(h_full));
postNoiseEnd = length(h_full);
if postNoiseEnd > postNoiseStart
    postNoise = h_full(postNoiseStart:postNoiseEnd);
else
    postNoise = [];
end

% 合并噪声估计
allNoise = [preNoise(:); postNoise(:)];
if ~isempty(allNoise)
    % 使用绝对中位差(MAD)估计噪声，对异常值更鲁棒
    noiseEst = 1.4826 * median(abs(allNoise - median(allNoise)));
    noiseEnergy = noiseEst^2 * length(signalSegment); % 标量到向量
else
    noiseEnergy = 1e-12;
end

% SNR计算
if noiseEnergy > 0 && signalEnergy > 0
    snr = 10*log10(signalEnergy / noiseEnergy);
else
    snr = -10; % 无效值
end

% ============== 可靠性判定 ==============
% 峰值突出度
if noiseEst > 0
    peakProminence = peakValLocal / noiseEst;
else
    peakProminence = 0;
end
peakHighEnough = (peakProminence > 5); % 降低阈值提高灵敏度

% 能量集中度
totalEnergy = sum(h_abs.^2);
energyConcentration = signalEnergy / (totalEnergy + 1e-12);
energyConcentrated = (energyConcentration > 0.03); % 降低阈值

% 预回声检查
preEchoEnd = max(1, peakIdx - 10); % 峰值前10个样本
if preEchoEnd > 1
    preEchoEnergy = sum(h_abs(1:preEchoEnd).^2);
    preEchoRatio = preEchoEnergy / (totalEnergy + 1e-12);
    lowPreEcho = (preEchoRatio < 0.3); % 放宽阈值
else
    lowPreEcho = true;
end

% 综合可靠性
peakReliable = peakHighEnough && energyConcentrated && lowPreEcho;

% ============== IR裁剪与后处理 ==============
% 确定IR起始点
irStart = max(1, peakIdx - cfg.deconvPreDelayKeep);
irEnd = min(irStart + irMaxLen - 1, length(h_full));

if irEnd >= irStart
    irRaw = h_full(irStart:irEnd);
else
    % 异常情况处理
    irRaw = zeros(min(irMaxLen, length(h_full)), 1);
    irStart = 1;
    irEnd = length(irRaw);
end

% 补零或截断到目标长度
if length(irRaw) < irMaxLen
    irOut = [irRaw; zeros(irMaxLen - length(irRaw), 1)];
else
    irOut = irRaw(1:irMaxLen);
end

% 调整峰值位置
peakPos = peakIdx - irStart + 1;
if peakPos < 1
    peakPos = 1;
elseif peakPos > irMaxLen
    peakPos = irMaxLen;
end

% ============== 返回结果 ==============
result = struct(...
    'ir', irOut, ...
    'peakPos', peakPos, ...
    'snr', snr, ...
    'peakVal', peakValLocal, ...
    'noiseEst', noiseEst, ...
    'signalEnergy', signalEnergy, ...
    'noiseEnergy', noiseEnergy, ...
    'energyConcentration', energyConcentration, ...
    'preEchoRatio', preEchoRatio, ...
    'peakProminence', peakProminence, ...
    'reliable', peakReliable, ...
    'irStartIdx', irStart, ...
    'irEndIdx', irEnd, ...
    'peakIdx', peakIdx);
end