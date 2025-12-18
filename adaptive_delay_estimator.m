function [delayD_new, st_delay, info] = adaptive_delay_estimator(st_delay, refFrame, errFrame, cfg)
% 自适应延迟估计器（基于 GCC-PHAT / 帧级互相关 + 平滑/约束）
% 
% 接口:
%   st_delay: 延迟估计状态（首次调用可传空 []，将被初始化）
%   refFrame: [N x nr] 参考通道帧（一般为 refUsed）
%   errFrame: [N x ne] 误差通道帧（一般为 errRaw）
%   cfg: 配置（必须包含 fs, minDelaySamples, maxDelaySamples, adaptiveDelayWindowSec,
%        adaptiveDelayUpdateFrames, adaptiveDelayConfidenceThreshold,
%        adaptiveDelaySmoothAlpha, adaptiveDelayMaxChangeSamples）
%
% 返回:
%   delayD_new: 建议的 delayD（样本）
%   st_delay: 更新后的状态结构
%   info: 估计信息结构，包含 estSamples, confidence, smoothed, changed (bool)
%
% 说明：
% - 采用短时 GCC-PHAT（FFT 权重归一化）在滑动窗口上估计参考与误差的相对时延。
% - 使用平滑滤波 (EMA) 与最大变化步长约束以避免突变。
% - 只在满足置信度阈值并且达到更新周期时采纳新的估计。
% - 对多通道输入默认做通道平均来获得单通道参考/误差用于估计。

% ------------- 初始化默认参数 -------------
fs = cfg.fs;
epsVal = 1e-12;

% ------------- 状态初始化 -------------
winLen = max(1, round(cfg.adaptiveDelayWindowSec * fs)); % 滑动窗口长度（样本）
if isempty(st_delay)
    st_delay = struct();
    st_delay.fs = fs;
    st_delay.winLen = winLen;
    st_delay.bufRef = zeros(winLen,1);
    st_delay.bufErr = zeros(winLen,1);
    st_delay.counter = 0;
    st_delay.smoothed = round(getfieldOrDefault(cfg, 'initialDelaySamples', cfg.minDelaySamples));
    st_delay.lastAccepted = st_delay.smoothed;
end

% ------------- 将多通道变为单通道（默认均值/或可定制）-------------
% 如果是多通道，先平均（你可以改为选主通道或波束成形后传入）
refMono = mean(refFrame, 2);
errMono = mean(errFrame, 2);

% ------------- 更新循环缓冲（简单移位拼接，保证长度恒定）-------------
N = length(refMono);
if N >= st_delay.winLen
    % 若帧长>=窗口长度，直接取最后 winLen 样本
    st_delay.bufRef = refMono(end-st_delay.winLen+1:end);
    st_delay.bufErr = errMono(end-st_delay.winLen+1:end);
else
    % 左移并追加
    st_delay.bufRef = [st_delay.bufRef(N+1:end); refMono];
    st_delay.bufErr = [st_delay.bufErr(N+1:end); errMono];
end

% ------------- 更新计数器，按周期触发估计 -------------
st_delay.counter = st_delay.counter + 1;
delayD_new = st_delay.smoothed;
info = struct('estSamples', [], 'confidence', 0, 'smoothed', st_delay.smoothed, 'changed', false);

if mod(st_delay.counter, cfg.adaptiveDelayUpdateFrames) ~= 0
    return; % 未到更新周期
end

% ------------- 做 GCC-PHAT（频域归一化互相关） -------------
x = st_delay.bufRef(:);
y = st_delay.bufErr(:);
L = length(x);

% FFT 长度（取 2^p >= 2*L）
nfft = 2^nextpow2(2*L);
X = fft(x, nfft);
Y = fft(y, nfft);
R = X .* conj(Y);
% PHAT 归一化
R = R ./ (abs(R) + epsVal);
r = real(ifft(R)); % circular cross-correlation
% 重新排列成 -nfft/2..nfft/2-1 的 lag 序列
r = fftshift(r);
lags = (-nfft/2):(nfft/2-1);

% ------------- 限制搜索范围（只考虑合理的延迟范围）-------------
maxLag = min(cfg.maxDelaySamples, round(cfg.adaptiveDelayWindowSec*fs)); % 不超过窗口长度
validIdx = find(abs(lags) <= maxLag);
[peakVal, idxRel] = max(abs(r(validIdx)));
idxGlobal = validIdx(idxRel);
lagAtPeak = lags(idxGlobal); % lag in samples (positive or negative)

% ------------- 置信度度量（峰值与旁瓣比 / 能量比）-------------
medianSide = median(abs(r(validIdx)));
if medianSide < epsVal, medianSide = epsVal; end
confidence = peakVal / medianSide; % ratio
% ------------- 得到估计（符号约定说明见下文）-------------
% 解释: 我们使用 r = IFFT( PHAT(X * conj(Y)) )，若 r 在正 lag = k 处最大，
% 表示 x[n - k] 与 y[n] 相关较强 => x leads y by k samples.
% 因此，若要对齐 ref -> err，估计的 delay (ref -> err) = k.
% 我们设置:
estimatedLag = lagAtPeak;
estimatedDelay = estimatedLag; % positive => ref leads err by estimatedDelay samples

% ------------- 平滑与约束 -------------
% EMA 平滑
smoothAlpha = cfg.adaptiveDelaySmoothAlpha;
smoothed = round(smoothAlpha * estimatedDelay + (1 - smoothAlpha) * st_delay.smoothed);

% 限制变化幅度
maxChange = cfg.adaptiveDelayMaxChangeSamples;
delta = smoothed - st_delay.lastAccepted;
if abs(delta) > maxChange
    smoothed = st_delay.lastAccepted + sign(delta) * maxChange;
end

% 忽略微小变化
if abs(smoothed - st_delay.lastAccepted) < cfg.adaptiveDelayMinChangeSamples
    info.changed = false;
else
    info.changed = true;
end

% 置信度门控
if confidence < cfg.adaptiveDelayConfidenceThreshold
    % 置信度不足，不接受本次估计（但仍更新 smoothed 状态用于下次）
    accepted = false;
else
    accepted = true;
end

% ------------- 最终采纳并约束在 min/max 范围 -------------
smoothed_clamped = min(max(smoothed, cfg.minDelaySamples), cfg.maxDelaySamples);
if accepted
    delayD_new = smoothed_clamped;
    st_delay.lastAccepted = delayD_new;
else
    delayD_new = st_delay.lastAccepted; % 保持上一次被接受的值
end

% 更新状态
st_delay.smoothed = smoothed; 
info.estSamples = estimatedDelay;
info.confidence = confidence;
info.smoothed = st_delay.smoothed;
info.accepted = accepted;
info.peakVal = peakVal;
info.lagAtPeak = lagAtPeak;

end

%% 小工具
function val = getfieldOrDefault(s, field, default)
    if isstruct(s) && isfield(s, field), val = s.(field); else val = default; end
end