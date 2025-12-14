function check_secondary_path(filePath)
% CHECK_SECONDARY_PATH 检查次级路径文件的可靠性
%   filePath: secondary_path.mat 文件的路径

if nargin < 1
    filePath = 'secondary_path/secondary_path.mat';
end

fprintf('🔍 正在加载次级路径文件: %s\n', filePath);
sec = load(filePath);

% === 1. 基础信息 ===
fprintf('\n--- 📋 基础信息 ---\n');
fprintf('采样率: %d Hz\n', sec.secondary.fs);
fprintf('扬声器数量: %d\n', sec.secondary.numSpeakers);
fprintf('麦克风数量: %d\n', sec.secondary.numMics);
fprintf('IR长度: %d 样本\n', sec.secondary.irLength);
fprintf('录制时间: %s UTC\n', sec.secondary.timestampUtc);

% === 2. 推荐延迟向量分析 (关键!) ===
delayVec = sec.secondary.delayEstimateSamples;
fprintf('\n--- ⏱️ 推荐延迟向量分析 ---\n');
fprintf('推荐延迟 (样本): [%s]\n', num2str(delayVec));
fprintf('推荐延迟 (ms):    [%s]\n', num2str(delayVec / sec.secondary.fs * 1000, '%.2f'));

% 检查是否过于整齐
if all(diff(delayVec) == 0)
    fprintf('⚠️  警告: 所有扬声器的推荐延迟完全相同! 这通常不正常。\n');
else
    fprintf('✅ 各扬声器延迟不同，符合物理预期。\n');
end

% === 3. 原始IR峰值位置检查 (最真实!) ===
ir = sec.secondary.impulseResponses;
[~, numMics, numSpks] = size(ir);
fprintf('\n--- 🗻 原始IR峰值位置 (绝对样本) ---\n');
allPeaks = zeros(numMics, numSpks);
for s = 1:numSpks
    for m = 1:numMics
        [~, pk] = max(abs(ir(:, m, s)));
        allPeaks(m, s) = pk;
        fprintf('Spk%d → Mic%d: %d 样本 (%.2f ms)\n', s, m, pk, pk/sec.secondary.fs*1000);
    end
end

% 计算每扬声器的延迟范围
fprintf('\n--- 📊 每扬声器延迟统计 ---\n');
for s = 1:numSpks
    peaks_s = allPeaks(:, s);
    range_s = max(peaks_s) - min(peaks_s);
    fprintf('Spk%d: 最小=%d, 最大=%d, 范围=%d 样本 (%.2f ms)\n', ...
        s, min(peaks_s), max(peaks_s), range_s, range_s/sec.secondary.fs*1000);
    if range_s > 500 % @48kHz, ~10ms 差异很大
        fprintf('    ⚠️  警告: Spk%d 到不同麦克风的延迟差异过大!\n', s);
    end
end

% === 4. 元数据分析 ===
fprintf('\n--- 📈 元数据分析 ---\n');
meta = sec.secondary.meta;
for s = 1:numSpks
    sMeta = meta.perSpeaker{s};
    fprintf('Spk%d:\n', s);
    fprintf('  - 可用性 (usable): %d\n', sMeta.usable);
    fprintf('  - 可靠峰值数: %d\n', sMeta.reliableCount);
    fprintf('  - 可靠峰值中位数: %d\n', sMeta.reliableMedian);
    fprintf('  - IQR: %s\n', num2str(sMeta.reliableIQR));
    fprintf('  - 尾部噪声 RMS: %.3e\n', sMeta.tailRMS);
    
    if ~sMeta.usable
        fprintf('    ❌ 警告: Spk%d 被标记为不可用!\n', s);
    end
end

% === 5. 可视化 (可选) ===
if ishandle(0) % 简单判断是否在图形环境
    choice = input('是否显示IR波形图? (y/n): ', 's');
    if strcmpi(choice, 'y')
        visualize_ir(ir, sec.secondary.fs, numSpks, numMics);
    end
end

fprintf('\n✅ 检查完成!\n');

end

function visualize_ir(ir, fs, numSpks, numMics)
    t_ms = (0:size(ir,1)-1) / fs * 1000;
    figure('Name', 'Secondary Path Impulse Responses', 'NumberTitle', 'off');
    for s = 1:numSpks
        for m = 1:numMics
            subplot(numSpks, numMics, (s-1)*numMics + m);
            plot(t_ms, ir(:, m, s), 'LineWidth', 1.2);
            grid on;
            xlabel('Time (ms)');
            ylabel('Amplitude');
            title(sprintf('Spk%d → Mic%d', s, m));
            % 标记峰值
            [~, pk_idx] = max(abs(ir(:, m, s)));
            hold on;
            plot(t_ms(pk_idx), ir(pk_idx, m, s), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            hold off;
        end
    end
    sgtitle('Secondary Path Impulse Responses with Peak Markers');
end