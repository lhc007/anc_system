function generate_quality_report(secondary, meta, cfg)
% 生成质量报告
% 包括统计数据和可视化

fprintf('\n[report] 生成质量报告...\n');

reportFile = strrep(cfg.secondaryPathFile, '.mat', '_report.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '次级路径测量质量报告\n');
fprintf(fid, '=======================\n\n');
fprintf(fid, '生成时间: %s\n', datestr(now));
fprintf(fid, '采样率: %d Hz\n', secondary.fs);
fprintf(fid, 'IR长度: %d 样本\n', secondary.irLength);
fprintf(fid, '扬声器数量: %d\n', secondary.numSpeakers);
fprintf(fid, '误差麦克风数量: %d\n', secondary.numMics);
fprintf(fid, '重复次数: %d\n\n', cfg.repetitions);

fprintf(fid, '各扬声器质量统计:\n');
fprintf(fid, '%-10s %-10s %-10s %-10s %-10s\n', ...
    '扬声器', '可用', 'SNR(dB)', '稳定性', '峰值位置');
fprintf(fid, '%s\n', repmat('-', 1, 60));

for spk = 1:cfg.numSpeakers
    if isfield(meta, 'perSpeaker') && length(meta.perSpeaker) >= spk
        spkMeta = meta.perSpeaker{spk};
        
        if isfield(spkMeta, 'qualityMetrics')
            metrics = spkMeta.qualityMetrics;
            peakPos = mean(spkMeta.peakPos(:));
            
            fprintf(fid, '%-10d %-10s %-10.1f %-10.2f %-10.0f\n', ...
                spk, ternary(metrics.usable, '是', '否'), ...
                metrics.medianSNR, metrics.stability, peakPos);
        else
            fprintf(fid, '%-10d %-10s %-10s %-10s %-10s\n', ...
                spk, '未知', 'N/A', 'N/A', 'N/A');
        end
    end
end

fprintf(fid, '\n质量指标说明:\n');
fprintf(fid, '1. SNR: 信噪比，应大于%.1f dB\n', cfg.snrThresholdDB);
fprintf(fid, '2. 稳定性: 峰值位置重复性，1.0为最佳\n');
fprintf(fid, '3. 峰值位置: IR主峰位置（样本）\n');

fclose(fid);
fprintf('    报告已保存: %s\n', reportFile);

% 生成可视化（如果MATLAB有图形界面）
if usejava('desktop')
    try
        visualize_results(secondary, meta, cfg);
        fprintf('    可视化图表已生成\n');
    catch ME
        fprintf('    可视化失败: %s\n', ME.message);
    end
end
end

function visualize_results(secondary, meta, cfg)
% 结果可视化
% 生成IR波形图和频响图

figure('Position', [100, 100, 1200, 800]);

% 子图1: 时域IR
subplot(2, 2, 1);
hold on;
colors = lines(cfg.numSpeakers);

for spk = 1:cfg.numSpeakers
    if meta.perSpeaker{spk}.usable
        ir = squeeze(secondary.impulseResponses(:, 1, spk));
        time = (0:length(ir)-1)/secondary.fs * 1000; % ms
        plot(time, ir, 'Color', colors(spk,:), 'LineWidth', 1.5);
    end
end

xlabel('时间 (ms)');
ylabel('幅度');
title('脉冲响应（时域）');
grid on;
legend(arrayfun(@(x) sprintf('Spk%d', x), 1:cfg.numSpeakers, 'UniformOutput', false));

% 子图2: 频响
subplot(2, 2, 2);
hold on;

for spk = 1:cfg.numSpeakers
    if meta.perSpeaker{spk}.usable
        ir = squeeze(secondary.impulseResponses(:, 1, spk));
        [H, F] = freqz(ir, 1, 1024, secondary.fs);
        semilogx(F, 20*log10(abs(H)+1e-12), 'Color', colors(spk,:), 'LineWidth', 1.5);
    end
end

xlabel('频率 (Hz)');
ylabel('幅度 (dB)');
title('频率响应');
grid on;
xlim([20, secondary.fs/2]);

% 子图3: SNR统计
subplot(2, 2, 3);
snrValues = [];
spkLabels = {};

for spk = 1:cfg.numSpeakers
    if isfield(meta.perSpeaker{spk}, 'qualityMetrics')
        snrValues(end+1) = meta.perSpeaker{spk}.qualityMetrics.medianSNR;
        spkLabels{end+1} = sprintf('Spk%d', spk);
    end
end

bar(snrValues);
hold on;
plot(xlim, [cfg.snrThresholdDB, cfg.snrThresholdDB], 'r--', 'LineWidth', 2);
ylabel('SNR (dB)');
title('各扬声器SNR');
set(gca, 'XTickLabel', spkLabels);
grid on;

% 子图4: 可用性统计
subplot(2, 2, 4);
usableCount = sum(cellfun(@(x) x.usable, meta.perSpeaker));
pie([usableCount, cfg.numSpeakers - usableCount], ...
    {'可用', '不可用'});
title(sprintf('可用性统计 (%d/%d)', usableCount, cfg.numSpeakers));

% 保存图像
saveDir = fileparts(cfg.secondaryPathFile);
saveas(gcf, fullfile(saveDir, 'secondary_path_visualization.png'));
end