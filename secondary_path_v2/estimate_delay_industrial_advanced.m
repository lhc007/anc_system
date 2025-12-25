function [delay, quality] = estimate_delay_industrial_advanced(recorded, sweep, cfg)
% 络延迟估计
% 质量评估

% === 强制对齐输入信号长度 ===
recorded = recorded(:); % 确保为列向量
sweep = sweep(:);

L_sweep = length(sweep);
L_recorded = length(recorded);

% 将较长的信号截断，或将较短的信号补零，使其长度一致
if L_recorded > L_sweep
    recorded = recorded(1:L_sweep);
elseif L_recorded < L_sweep
    recorded = [recorded; zeros(L_sweep - L_recorded, 1)];
end
% 现在 recorded 和 sweep 长度相同，可以安全使用 'coeff' 选项

% 计算能量包络
energy_rec = recorded.^2;
energy_sweep = sweep.^2;
fs = cfg.fs;
% 简单平滑（移动平均）
window_len = round(0.01 * fs);  % 10ms窗口
window = ones(window_len, 1)/window_len;
env_rec = conv(energy_rec, window, 'same');
env_sweep = conv(energy_sweep, window, 'same');

% 降采样
decim_factor = max(1, floor(fs / 2000));  % 降至2kHz
env_rec_ds = downsample(env_rec, decim_factor);
env_sweep_ds = downsample(env_sweep, decim_factor);
fs_ds = fs / decim_factor;

% 互相关
max_lag = ceil(0.5 * fs_ds);
if length(env_rec_ds) > 2*max_lag && length(env_sweep_ds) > 2*max_lag
    [C, lags] = xcorr(env_rec_ds, env_sweep_ds, max_lag, 'coeff');
else
    [C, lags] = xcorr(env_rec_ds, env_sweep_ds, 'coeff');
end

% 寻找主峰
[peak_val, peak_idx] = max(C);
delay_ds = lags(peak_idx);
delay = delay_ds * decim_factor;

% 改进的质量评估
% 1. 峰值高度
height_score = min(1.0, peak_val * 1.5);

% 2. 峰值锐度（寻找次高峰）
C_sorted = sort(C, 'descend');
if length(C_sorted) >= 2
    prominence = (C_sorted(1) - C_sorted(2)) / (C_sorted(1) + eps);
    prominence_score = min(1.0, prominence * 3);
else
    prominence_score = 0.5;
end

% 3. 搜索范围内的能量集中度
peak_range = max(1, peak_idx-5):min(length(C), peak_idx+5);
energy_in_peak = sum(C(peak_range).^2);
total_energy = sum(C.^2);
energy_score = energy_in_peak / (total_energy + eps);

% 综合质量评分
quality = 0.4*height_score + 0.4*prominence_score + 0.2*energy_score;
quality = max(0.1, min(1.0, quality));

% 确保延迟为正
delay = max(cfg.minPhysDelaySamples, ...
    min(cfg.maxPhysDelaySamples, delay));
end