%% check_acoustic_signal_strength.m
% 用于评估次级路径测量中实际接收到的声学信号强度和 SNR

clear; clc;

% === 配置 ===
noise_rms_est = 1.9e-4;  % 你实测的麦克风本底噪声 RMS（来自 check_noise_floor.m）
window_half_len = 200;   % 主峰前后取多少样本计算信号能量
secondary_file = 'secondary_path/secondary_path.mat';

% === 加载数据 ===
if ~exist(secondary_file, 'file')
    error('未找到次级路径文件: %s', secondary_file);
end
load(secondary_file, 'secondary');

numSpk = secondary.numSpeakers;
numMic = secondary.numMics;
fs = secondary.fs;

fprintf('采样率: %d Hz\n', fs);
fprintf('扬声器数量: %d, 麦克风数量: %d\n', numSpk, numMic);

% === 分析每个 Spk → Mic 通道 ===
figure('Name', '冲激响应与信号强度分析', 'Position', [100, 100, 1200, 800]);
t = (0:size(secondary.impulseResponses,1)-1) / fs * 1000; % 时间轴 (ms)

for spk = 1:numSpk
    for mic = 1:numMic
        ir = secondary.impulseResponses(:, mic, spk);
        
        % 找主峰位置（最大绝对值）
        [~, peakIdx] = max(abs(ir));
        
        % 定义信号窗口（主峰附近）
        win_start = max(1, peakIdx - window_half_len);
        win_end   = min(length(ir), peakIdx + window_half_len);
        signal_segment = ir(win_start:win_end);
        
        % 计算信号 RMS
        signal_rms = rms(signal_segment);
        
        % 计算 SNR (dB)
        snr_db = 20 * log10(signal_rms / noise_rms_est);
        
        % 显示结果
        subplot(numSpk, numMic, (spk-1)*numMic + mic);
        plot(t, ir, 'b', 'LineWidth', 1);
        hold on;
        % 标出主峰
        plot(t(peakIdx), ir(peakIdx), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        % 标出信号窗口
        yl = ylim;
        fill([t(win_start) t(win_end) t(win_end) t(win_start)], ...
             [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor', 'none');
        hold off;
        
        title(sprintf('Spk%d → Mic%d\nSNR=%.1f dB', spk, mic, snr_db), 'FontSize', 9);
        xlabel('时间 (ms)'); ylabel('幅度');
        grid on;
        xlim([0, min(200, t(end))]); % 只看前 200ms
        
        % 打印到命令行
        fprintf('Spk%d → Mic%d: 信号RMS=%.2e, SNR=%.1f dB\n', ...
                spk, mic, signal_rms, snr_db);
    end
end

sgtitle('次级路径冲激响应与实测 SNR 分析', 'FontSize', 12);

% === 总体判断 ===
all_snrs = [];
for spk = 1:numSpk
    for mic = 1:numMic
        ir = secondary.impulseResponses(:, mic, spk);
        [~, peakIdx] = max(abs(ir));
        win = max(1, peakIdx-200):min(length(ir), peakIdx+200);
        sig_rms = rms(ir(win));
        snr = 20*log10(sig_rms / noise_rms_est);
        all_snrs(end+1) = snr;
    end
end

mean_snr = mean(all_snrs);
fprintf('\n>>> 平均 SNR = %.1f dB\n', mean_snr);
if mean_snr > 10
    fprintf('✅ 声学信号强度足够，次级路径可能可用。\n');
else
    fprintf('❌ 声学信号太弱！请提高扬声器音量或改善声学耦合。\n');
end