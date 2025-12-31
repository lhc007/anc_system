function [sweepCore_scaled, sweepSig, info] = generate_excitation_signal(cfg)
% generate_excitation_signal - 生成激励信号（扫频信号）
% 构建播放信号 [前静音][扫频][后静音]
% 包含完整的ESS生成和信号构建流程

fs = cfg.fs;

% 生成原始ESS
[sweepCore, ~, ~, info] = generate_ess_core(cfg);

% 构建播放信号 [前静音][扫频][后静音]
pre_silence_samples = round(cfg.padLeading * fs);
post_silence_samples = round(cfg.padTrailing * fs);

% 构建信号：[pre_silence][sweepCore][post_silence]
sweepCore_scaled = sweepCore(:);
if max(abs(sweepCore_scaled)) > 0
    sweepCore_scaled = (sweepCore_scaled / max(abs(sweepCore_scaled))) * cfg.amplitude;
end

sweepSig = [ ...
    zeros(pre_silence_samples, 1); ...
    sweepCore_scaled; ...
    zeros(post_silence_samples, 1) ...
];

% 确保列向量
sweepSig = sweepSig(:);    % 播放信号
end

function [sweepCore, sweepFull, invFilter, info] = generate_ess_core(cfg)
% generate_ess_core - 生成 Exponential Sine Sweep (ESS) 核心信号（不含任何静音段）
%
% 输出:
%   sweepCore  : 纯扫频信号（未归一化，长度 = round(sweepDuration * fs)）
%   sweepFull  : 完整激励信号（含前后静音，长度 = N + padLead + padTrail）
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

% 生成核心扫频
t = (0:N-1)' / fs;
% 指数扫频相位
phase = 2 * pi * f1 * T / log(f2/f1) * (exp(t * log(f2/f1) / T) - 1);
sweepCore = sin(phase);

% 构造反卷积滤波器（幅度加权倒谱）
% w(t) = exp(-t * log(f2/f1) / T)
weight = exp(-t * log(f2/f1) / T);
invFilter = flipud(weight .* sweepCore);  % 时间反转 + 加权


% === 生成完整激励信号（含前后静音）===
padLeadSamples = round(cfg.padLeading * fs);
padTrailSamples = round(cfg.padTrailing * fs);

% 拼接完整信号：[前静音][核心扫频][后静音]
sweepFull = [zeros(padLeadSamples, 1); 
             sweepCore; 
             zeros(padTrailSamples, 1)];

% 可选：预加重/低频增强（此处不处理，交由上层决定）

% 返回元信息
info.N = N;
info.f1 = f1;
info.f2 = f2;
info.fs = fs;
info.sweepDuration = T;
info.padLeading = cfg.padLeading;
info.padTrailing = cfg.padTrailing;
info.totalLength = length(sweepFull);

end