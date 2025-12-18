%% check_acoustic_signal_strength.m
% 用于评估次级路径测量中实际接收到的声学信号强度和 SNR
% ✅ 修复：正确显示物理麦克风通道号（如 Mic5, Mic6）

clear; clc;
cfg = anc_config();
% === 配置 ===
noise_rms_est = estimate_noise_floor('Duration', 2.0, 'UseCache', false, 'Verbose', true);
% 推荐：固定时间窗 25 ms
analysis_window_ms = 25;
window_half_len = round((cfg.fs * analysis_window_ms / 1000) / 2);
secondary_file = 'secondary_path/secondary_path.mat';

% === 加载数据 ===
if ~exist(secondary_file, 'file')
    error('未找到次级路径文件: %s', secondary_file);
end
load(secondary_file, 'secondary');

numSpk = secondary.numSpeakers;
numMic = secondary.numMics;
fs = secondary.fs;

% === 获取物理麦克风通道号 ===
if isfield(secondary, 'errorMicPhysicalChannels')
    micPhysical = secondary.errorMicPhysicalChannels(:)';  % 行向量，如 [5, 6]
    if length(micPhysical) ~= numMic
        warning('物理通道数与保存的麦克风数不一致，使用默认编号 1:%d', numMic);
        micPhysical = 1:numMic;
    end
else
    % 尝试从配置回退
    if isfield(secondary, 'measureConfig') && ...
       isfield(secondary.measureConfig, 'micChannels') && ...
       isfield(secondary.measureConfig.micChannels, 'error')
        allErr = secondary.measureConfig.micChannels.error;
        micPhysical = allErr(1:numMic);
    else
        micPhysical = 1:numMic;
        fprintf('⚠️ 未找到物理通道信息，假设麦克风为: %s\n', mat2str(micPhysical));
    end
end

fprintf('采样率: %d Hz\n', fs);
fprintf('扬声器数量: %d, 麦克风数量: %d\n', numSpk, numMic);
fprintf('物理误差麦克风通道: %s\n', mat2str(micPhysical));

% === 分析每个 Spk → Mic 通道 ===
figure('Name', '冲激响应与信号强度分析', 'Position', [100, 100, 1200, 800]);
t = (0:size(secondary.impulseResponses,1)-1) / fs * 1000; % 时间轴 (ms)

all_snrs = [];

for spk = 1:numSpk
    for mi = 1:numMic
        ir = secondary.impulseResponses(:, mi, spk);  % 逻辑索引 mi = 1,2
        
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
        all_snrs(end+1) = snr_db;
        
        % 显示结果（使用物理通道号）
        subplot(numSpk, numMic, (spk-1)*numMic + mi);
        plot(t, ir, 'b', 'LineWidth', 1);
        hold on;
        plot(t(peakIdx), ir(peakIdx), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        yl = ylim;
        fill([t(win_start) t(win_end) t(win_end) t(win_start)], ...
             [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor', 'none');
        hold off;
        
        title(sprintf('Spk%d → Mic%d\nSNR=%.1f dB', spk, micPhysical(mi), snr_db), 'FontSize', 9);
        xlabel('时间 (ms)'); ylabel('幅度');
        grid on;
        xlim([0, min(200, t(end))]); % 只看前 200ms
        
        % 打印到命令行（使用物理通道号）
        fprintf('Spk%d → Mic%d: 信号RMS=%.2e, SNR=%.1f dB\n', ...
                spk, micPhysical(mi), signal_rms, snr_db);
    end
end

sgtitle('次级路径冲激响应与实测 SNR 分析（物理麦克风通道）', 'FontSize', 12);

% === 总体判断 ===
mean_snr = mean(all_snrs);
fprintf('\n>>> 平均 SNR = %.1f dB\n', mean_snr);
if mean_snr > 10
    fprintf('✅ 声学信号强度足够，次级路径可能可用。\n');
elseif mean_snr > 5
    fprintf('⚠️ 信号较弱，ANC性能可能受限。建议提高音量或靠近扬声器。\n');
else
    fprintf('❌ 声学信号太弱！请提高扬声器音量或改善声学耦合。\n');
end