function snr_db = estimate_snr(y_rec, sweepSig, fs)
% ESTIMATE_SNR: 基于激励信号尾部静音段估计 SNR
%
% 输入:
%   y_rec      : [N x M] 录制的麦克风信号（包含激励响应 + 后续静音）
%   sweepSig   : [L x 1] 激励信号（必须与播放的信号完全一致，含 padTrailing 静音）
%   fs         : 标量，采样率 (Hz)
%   cfg        : 配置结构体，需包含 cfg.padTrailing（后静音时长，秒）
%
% 输出:
%   snr_db     : [1 x M] 各通道 SNR（dB），若信号太短则为 -Inf
%
% 假设: sweepSig 末尾有足够静音（由 cfg.padTrailing 保证）

[N, N_ch] = size(y_rec);
L_sweep = length(sweepSig);

% 找到 sweepSig 中最后一个非零样本（即激励结束点）
last_sig_idx = find(sweepSig ~= 0, 1, 'last');
if isempty(last_sig_idx)
    last_sig_idx = L_sweep;
end

% 噪声估计区域：从激励结束后开始，取 0.1 秒
noise_start_in_sweep = last_sig_idx + 1;
noise_duration_samples = round(0.1 * fs); % 100ms 噪声

snr_db = zeros(1, N_ch);

for ch = 1:N_ch
    y_ch = y_rec(:, ch);
    
    % 确保录制信号足够长以包含噪声段
    if length(y_ch) < L_sweep + noise_duration_samples
        warning('Channel %d: recorded signal too short for reliable SNR.', ch);
        snr_db(ch) = -Inf;
        continue;
    end
    
    % 噪声段：对应于 sweepSig 之后的静音部分
    noise_seg = y_ch(L_sweep + (1:noise_duration_samples));
    
    % 信号段：整个 sweep 响应（可用 ir_est 能量，但此处用录制信号近似）
    signal_seg = y_ch(1:L_sweep);
    
    % 计算功率
    noise_power = mean(noise_seg.^2);
    signal_power = mean(signal_seg.^2);
    
    if noise_power <= 0
        snr_db(ch) = Inf;
    else
        snr_db(ch) = 10 * log10(signal_power / noise_power);
    end
end
end