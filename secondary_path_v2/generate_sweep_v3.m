function [sweepCore, invFilter, info] = generate_sweep_v3(cfg)
% GENERATE_SWEEP - 生成 Exponential Sine Sweep (ESS) 核心信号（不含任何静音段）
%
% 输出:
%   sweepCore  : 纯扫频信号（未归一化，长度 = round(sweepDuration * fs)）
%   invFilter  : 对应的反卷积滤波器（用于后续 IR 提取）
%   info       : 结构体，包含 f1, f2, N, fs 等元信息

fs = cfg.fs;
f1 = cfg.sweepF1;
f2 = cfg.sweepF2;
T = cfg.sweepDuration;               % 扫频持续时间（秒）

% 计算样本数（仅核心部分）
N = round(T * fs);
if N <= 0
    error('无效的 sweepDuration 或采样率');
end

t = (0:N-1)' / fs;

% 指数扫频相位
phase = 2 * pi * f1 * T / log(f2/f1) * (exp(t * log(f2/f1) / T) - 1);

% 生成核心扫频
sweepCore = sin(phase);

% 构造反卷积滤波器（幅度加权倒谱）
% w(t) = exp(-t * log(f2/f1) / T)
weight = exp(-t * log(f2/f1) / T);
invFilter = flipud(weight .* sweepCore);  % 时间反转 + 加权

% 可选：预加重/低频增强（此处不处理，交由上层决定）

% 返回元信息
info.N = N;
info.f1 = f1;
info.f2 = f2;
info.fs = fs;
info.sweepDuration = T;

end