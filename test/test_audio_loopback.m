% test_audio_loopback_optimized.m
% 优化版音频链路测试（无Audio Toolbox依赖）

clear; clc; close all;

fprintf('=== 音频链路测试（优化版）===\n\n');
cfg = anc_config();
% 基本配置
fs = cfg.fs;
duration = cfg.sweepDuration; % 延长测试时间到5秒
samples = round(fs * duration);
t = (0:samples-1)'/fs;

% 生成更丰富的测试信号
fprintf('1. 生成测试信号...\n');
test_signal = zeros(samples, 1);

% 第1部分：1kHz正弦波（0-2秒）
idx1 = t < 2;
test_signal(idx1) = 0.5 * sin(2*pi*1000*t(idx1));

% 第2部分：线性扫频500-2000Hz（2-4秒）
idx2 = (t >= 2) & (t < 4);
duration_sweep = 2;
t_sweep = t(idx2) - 2;
f_start = 500;
f_end = 2000;
sweep = chirp(t_sweep, f_start, duration_sweep, f_end, 'linear');
test_signal(idx2) = 0.5 * sweep;

% 第3部分：白噪声（4-5秒）
idx3 = t >= 4;
test_signal(idx3) = 0.2 * randn(sum(idx3), 1);

% 初始化硬件
hw = hardware_init_measure(cfg);

% 预热
fprintf('3. 预热音频设备...\n');
for i = 1:10
    hw.writer(zeros(cfg.timeFrameSamples, cfg.numSpeakers));
    hw.reader();
end

% 播放并录制
fprintf('4. 开始播放和录制（%d秒）...\n', duration);
fprintf('   测试顺序：1kHz正弦波 → 扫频信号 → 白噪声\n');
fprintf('   请确保扬声器正在播放测试音...\n');

numBlocks = ceil(samples / cfg.timeFrameSamples);
recorded = zeros(samples, cfg.micNumChannels);
ptr = 1;

for b = 1:numBlocks
    % 准备输出块
    nToCopy = min(cfg.timeFrameSamples, samples - ptr + 1);
    outBlock = zeros(cfg.timeFrameSamples, cfg.numSpeakers);
    
    if nToCopy > 0
        outBlock(1:nToCopy, 1) = test_signal(ptr:ptr+nToCopy-1);
    end
    
    % 播放并录制
    hw.writer(outBlock);
    micFrame = hw.reader();
    
    if ~isempty(micFrame) && size(micFrame,1) >= nToCopy
        if nToCopy > 0
            recorded(ptr:ptr+nToCopy-1, :) = micFrame(1:nToCopy, :);
        end
    end
    
    ptr = ptr + cfg.timeFrameSamples;
    
    % 显示进度
    if mod(b, 20) == 0
        progress = min(100, round(b/numBlocks*100));
        fprintf('   进度: %d%%\n', progress);
    end
end

% 释放硬件
hw.release();
fprintf('5. 音频播放和录制完成\n');

% 分析结果
fprintf('\n=== 分析结果 ===\n');

% 计算整体能量
overall_rms = rms(recorded(:));
fprintf('整体RMS值: %.4f\n', overall_rms);

% 分析每个通道
for ch = 1:cfg.micNumChannels
    fprintf('\n--- 通道 %d ---\n', ch);
    ch_signal = recorded(:, ch);
    
    % 基本统计
    rms_val = rms(ch_signal);
    max_val = max(abs(ch_signal));
    fprintf('RMS: %.4f, 峰值: %.4f\n', rms_val, max_val);
    
    % 分时段分析
    segment_duration = 1; % 每段1秒
    segment_samples = fs * segment_duration;
    num_segments = floor(samples / segment_samples);
    
    for seg = 1:min(5, num_segments) % 只分析前5秒
        start_idx = (seg-1)*segment_samples + 1;
        end_idx = seg*segment_samples;
        
        if end_idx <= samples
            segment = ch_signal(start_idx:end_idx);
            
            % 计算能量
            seg_energy = sum(segment.^2);
            
            % 如果是第一段（1kHz正弦波），计算1kHz能量占比
            if seg == 1
                % 计算频谱
                N = min(4096, length(segment));
                fft_data = fft(segment(1:N), N);
                freq = (0:N-1)' * fs / N;
                
                % 1kHz附近能量
                idx_1k = find(freq >= 950 & freq <= 1050);
                if ~isempty(idx_1k)
                    energy_1k = sum(abs(fft_data(idx_1k)).^2);
                    total_energy = sum(abs(fft_data).^2);
                    ratio_1k = energy_1k / total_energy;
                    
                    fprintf('  第1秒（1kHz正弦波）：能量=%.2e，1kHz占比=%.1f%%\n', ...
                        seg_energy, ratio_1k*100);
                end
            else
                fprintf('  第%d秒：能量=%.2e\n', seg, seg_energy);
            end
        end
    end
    
    % 判断通道状态
    if rms_val > 0.01
        fprintf('  状态: ✅ 强信号\n');
    elseif rms_val > 0.001
        fprintf('  状态: ⚠️  弱信号\n');
    else
        fprintf('  状态: ❌ 无信号\n');
    end
end

% 检查声道交叉
fprintf('\n=== 声道交叉分析 ===\n');
if cfg.micNumChannels >= 2
    % 计算通道1和通道2的相关性
    correlation = corrcoef(recorded(:,1), recorded(:,2));
    fprintf('通道1-2相关性: %.3f\n', correlation(1,2));
    
    if abs(correlation(1,2)) > 0.7
        fprintf('提示: 通道1和2高度相关，可能是相同的信号源\n');
    end
end

% 绘制结果
fprintf('\n=== 生成图表 ===\n');
figure('Position', [100, 100, 1400, 800]);

% 1. 时域对比图
subplot(3,3,1);
plot(t, test_signal, 'b', 'LineWidth', 1);
title('发送信号（扬声器）');
xlabel('时间 (s)');
ylabel('幅度');
xlim([0, duration]);
grid on;
legend('测试信号');

subplot(3,3,2);
plot(t, recorded(:,1:min(3,size(recorded,2))));
title('录制信号（前3个通道）');
xlabel('时间 (s)');
ylabel('幅度');
xlim([0, duration]);
grid on;
legend(arrayfun(@(x) sprintf('Ch%d', x), 1:min(3,size(recorded,2)), 'UniformOutput', false));

% 2. 频谱图（前3个通道）
subplot(3,3,3);
hold on;
colors = ['r', 'g', 'b', 'c', 'm', 'y'];
for ch = 1:min(3, cfg.micNumChannels)
    [pxx, f] = pwelch(recorded(:,ch), hamming(1024), 512, 1024, fs);
    plot(f, 10*log10(pxx), colors(ch), 'LineWidth', 1);
end
title('录制信号频谱');
xlabel('频率 (Hz)');
ylabel('功率 (dB)');
xlim([20, 5000]); % 限制到5kHz
grid on;
legend(arrayfun(@(x) sprintf('Ch%d', x), 1:min(3,cfg.micNumChannels), 'UniformOutput', false));

% 3. 时频分析（通道1）
subplot(3,3,4);
[s, f_spec, t_spec] = spectrogram(recorded(:,1), 256, 250, 256, fs, 'yaxis');
imagesc(t_spec, f_spec, 10*log10(abs(s)));
axis xy;
xlabel('时间 (s)');
ylabel('频率 (Hz)');
title('通道1时频图');
colorbar;
ylim([0, 5000]);
caxis([-80, 0]);

% 4. 发送和接收信号对比（前2秒，1kHz部分）
subplot(3,3,5);
idx_plot = t < 2;
plot(t(idx_plot), test_signal(idx_plot), 'b', 'LineWidth', 2);
hold on;
plot(t(idx_plot), recorded(idx_plot, 1), 'r--', 'LineWidth', 1);
title('1kHz正弦波对比');
xlabel('时间 (s)');
ylabel('幅度');
xlim([0, 0.02]); % 只看前20ms
grid on;
legend('发送', '接收 (Ch1)');

% 5. 各通道RMS值比较
subplot(3,3,6);
rms_values = zeros(cfg.micNumChannels, 1);
for ch = 1:cfg.micNumChannels
    rms_values(ch) = rms(recorded(:,ch));
end
bar(rms_values);
title('各通道RMS值');
xlabel('通道');
ylabel('RMS值');
grid on;
set(gca, 'XTick', 1:cfg.micNumChannels);

% 6. 相关性矩阵
subplot(3,3,7);
if cfg.micNumChannels > 1
    corr_matrix = corrcoef(recorded);
    imagesc(corr_matrix);
    colorbar;
    title('通道相关性矩阵');
    xlabel('通道');
    ylabel('通道');
    set(gca, 'XTick', 1:cfg.micNumChannels);
    set(gca, 'YTick', 1:cfg.micNumChannels);
    caxis([-1, 1]);
else
    text(0.5, 0.5, '需要多个通道', 'HorizontalAlignment', 'center');
    title('通道相关性矩阵');
end

% 7. 相位差分析（如果有多个通道）
subplot(3,3,8);
if cfg.micNumChannels >= 2
    % 计算通道1和其他通道的相位差
    phase_diffs = zeros(cfg.micNumChannels-1, 1);
    for ch = 2:cfg.micNumChannels
        [c, lags] = xcorr(recorded(:,1), recorded(:,ch));
        [~, max_idx] = max(abs(c));
        phase_diffs(ch-1) = lags(max_idx) / fs * 1000; % 转换为毫秒
    end
    
    bar(phase_diffs);
    title('相对于通道1的时延（毫秒）');
    xlabel('通道');
    ylabel('时延 (ms)');
    grid on;
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('Ch%d', x), 2:cfg.micNumChannels, 'UniformOutput', false));
else
    text(0.5, 0.5, '需要多个通道', 'HorizontalAlignment', 'center');
    title('相位差分析');
end

% 8. 信噪比估计
subplot(3,3,9);
snr_estimates = zeros(cfg.micNumChannels, 1);
for ch = 1:cfg.micNumChannels
    signal = recorded(:,ch);
    % 简单信噪比估计：假设最后0.5秒是噪声
    noise_segment = signal(end-round(0.5*fs)+1:end);
    noise_power = mean(noise_segment.^2);
    signal_power = mean(signal.^2);
    
    if noise_power > 0
        snr_estimates(ch) = 10*log10(signal_power/noise_power);
    else
        snr_estimates(ch) = Inf;
    end
end

bar(snr_estimates);
title('各通道信噪比估计');
xlabel('通道');
ylabel('SNR (dB)');
grid on;
set(gca, 'XTick', 1:cfg.micNumChannels);
ylim([0, max(50, max(snr_estimates(snr_estimates<Inf)))]);

fprintf('\n=== 测试完成 ===\n');

% 提供具体建议
fprintf('\n=== 诊断建议 ===\n');
if overall_rms < 0.001
    fprintf('❌ 严重问题：几乎没有录制到任何信号\n');
    fprintf('   可能原因：\n');
    fprintf('   1. 扬声器没有播放声音\n');
    fprintf('   2. 麦克风被静音或音量极低\n');
    fprintf('   3. 扬声器和麦克风之间距离太远或有障碍物\n');
    fprintf('\n   解决方案：\n');
    fprintf('   1. 检查Windows音量设置，确保扬声器2未被静音\n');
    fprintf('   2. 尝试用其他音频播放器测试扬声器\n');
    fprintf('   3. 将麦克风靠近扬声器测试\n');
    fprintf('   4. 检查Realtek Audio Console中的设置\n');
elseif overall_rms < 0.01
    fprintf('⚠️  警告：信号较弱\n');
    fprintf('   可能原因：\n');
    fprintf('   1. 扬声器音量太低\n');
    fprintf('   2. 麦克风灵敏度设置过低\n');
    fprintf('   3. 环境噪声干扰\n');
    fprintf('\n   解决方案：\n');
    fprintf('   1. 增加测试信号的幅度（当前为0.5）\n');
    fprintf('   2. 在Windows麦克风设置中提高增益\n');
    fprintf('   3. 在安静环境中测试\n');
else
    fprintf('✅ 信号强度正常\n');
    fprintf('   请检查频谱图是否显示1kHz和扫频信号\n');
end

% 保存数据以便进一步分析
save('audio_loopback_test.mat', 'test_signal', 'recorded', 't', 'fs', 'cfg');
fprintf('\n测试数据已保存为 audio_loopback_test.mat\n');