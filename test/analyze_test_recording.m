% analyze_test_recording_dual.m
% 分析双通道扫频响应（例如 CH5 & CH6），分别评估每个麦克风的可用性

clear; clc;
fs = 48000;

% 加载信号
[sweep, ~] = audioread('test_sweep_48k_5s.wav');
[recorded, fs_r] = audioread('recorded_ch5_ch6.wav');  % ← 双通道文件

if fs_r ~= fs
    error('采样率不匹配！请确保录制时使用 %d Hz', fs);
end

% 确保是双通道
if size(recorded, 2) ~= 2
    error('期望双通道录音（CH5 & CH6），但输入有 %d 通道', size(recorded, 2));
end

nChannels = 2;
channelNames = {'CH5', 'CH6'};

% 预分配结果
ir_all = {};
snr_db_all = zeros(1, nChannels);
peakIdx_all = zeros(1, nChannels);
maxAmp_all = zeros(1, nChannels);
usable_all = false(1, nChannels);

L = length(sweep);
Nfft = 2^nextpow2(2*L);
S = fft(sweep, Nfft);
InvFilter = conj(S) ./ (abs(S).^2 + 1e-12); % 正则化逆滤波器

% 对每个通道单独处理
% 在 analyze_test_recording_dual.m 中替换 for 循环内部：
for ch = 1:nChannels
    fprintf('\n处理 %s...\n', channelNames{ch});
    rec_sig = recorded(:, ch);
    
    [ir, snr_db, peakIdx, meta] = deconvolve_and_estimate_snr(...
        rec_sig, sweep, fs, ...
        'sigWinRadius', 200, ...
        'noiseWinStart', 2000, ...
        'noiseWinLength', 2000, ...
        'irMaxLenExtra', 10000);
    
    ir_all{ch} = ir;
    snr_db_all(ch) = snr_db;
    peakIdx_all(ch) = peakIdx;
    
    maxAmp_all(ch) = max(abs(rec_sig));
    usable_all(ch) = (snr_db > 10) && (peakIdx > 0.01*fs) && (maxAmp_all(ch) < 0.98);
end

% === 绘图 ===
figure('Position', [100, 100, 1200, 800]);

for ch = 1:nChannels
    ir = ir_all{ch};
    t_ir = (0:length(ir)-1) / fs;
    peakIdx = peakIdx_all(ch);
    
    subplot(2, 2, ch);  % 上排：原始 sweep（可选）
    plot((0:length(sweep)-1)/fs, sweep);
    title('原始扫频信号');
    xlabel('时间 (s)'); ylabel('幅度');
    grid on;
    
    subplot(2, 2, ch + 2);  % 下排：IR
    plot(t_ir, ir);
    hold on;
    plot(t_ir(peakIdx), ir(peakIdx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    title(sprintf('%s 脉冲响应 | 主峰=%.1f ms | SNR=%.1f dB', ...
        channelNames{ch}, 1000*t_ir(peakIdx), snr_db_all(ch)));
    xlabel('时间 (s)'); ylabel('幅度');
    grid on;
end

sgtitle('双通道扫频响应分析（CH5 & CH6）');

% === 打印评估结果 ===
fprintf('\n=== 麦克风可用性评估（双通道）===\n');
for ch = 1:nChannels
    t_peak_ms = 1000 * (peakIdx_all(ch) - 1) / fs;  % 样本索引从1开始，时间从0开始
    fprintf('\n--- %s ---\n', channelNames{ch});
    fprintf('主峰位置: %.1f ms (%d 样本)\n', t_peak_ms, peakIdx_all(ch));
    fprintf('估计 SNR: %.1f dB\n', snr_db_all(ch));
    fprintf('录制最大幅度: %.3f\n', maxAmp_all(ch));
    if usable_all(ch)
        fprintf('✅ %s 可用\n', channelNames{ch});
    else
        fprintf('❌ %s 不可用\n', channelNames{ch});
    end
end