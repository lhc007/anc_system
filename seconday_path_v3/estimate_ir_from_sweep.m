function [ir_est, ir_offsets] = estimate_ir_from_sweep(y_rec_crop, sweepSig, sweep_core, cfg)
% ESTIMATE_IR_FROM_SWEEP: 从裁剪后的录制信号中提取多通道次级路径冲激响应 (IR)
% 输入:
%   y_rec_crop : [L_crop × N_ch] 裁剪后的录制信号（已包含完整响应）
%   sweepSig   : [L_full × 1]    完整激励信号（含前后静音）
%   sweep_core : [T × 1]         核心扫频部分（用于验证，此处未使用）
%   cfg        : 配置结构体
% 输出:
%   ir_est     : [L_ir × N_ch]   估计的冲激响应
%   ir_offsets : [1 × N_ch]      IR 在 h_full 中的起始索引（局部坐标）

[M, N_ch] = size(y_rec_crop);
L_full = length(sweepSig);
L_ir = cfg.irMaxLen;
fs = cfg.fs;

% 预分配
ir_est = zeros(L_ir, N_ch);
ir_offsets = zeros(1, N_ch);
% coherence = zeros(1, N_ch);

% 构建逆滤波器（ESS 反卷积）
inv_filter = fliplr(conj(sweepSig)); 

for ch = 1:N_ch
    y_ch = y_rec_crop(:, ch);
    
    % === 1. 频域反卷积得到完整冲激响应 ===
    Nfft = 2^nextpow2(L_full + length(y_ch) - 1);
    H = fft(y_ch, Nfft) .* fft(inv_filter, Nfft);
    h_full = ifft(H, Nfft, 'symmetric');
    h_full = h_full(1:length(y_ch) + L_full - 1); % 精确线性卷积长度
    
    % === 2. 在合理范围内搜索主峰（10ms ～ 1.5s）===
    min_delay = round(0.01 * fs);   % 10ms
    max_delay = min(round(1.5 * fs), length(h_full));
    if max_delay <= min_delay
        max_delay = length(h_full);
    end
    
    [~, peak_local_idx] = max(abs(h_full(min_delay:max_delay)));
    peak_idx = min_delay + peak_local_idx - 1; % 主峰在 h_full 中的绝对索引
    
    % === 3. 提取固定长度 IR（保留 pre-delay）===
    ir_start = max(1, peak_idx - cfg.deconvPreDelayKeep);
    if ir_start + L_ir - 1 > length(h_full)
        ir_start = length(h_full) - L_ir + 1;
    end
    
    ir_est(:, ch) = h_full(ir_start : ir_start + L_ir - 1);
    ir_offsets(ch) = ir_start;
    
end
end