function [actual_start_idx, recorded_aligned] = align_sweep_start(recordedFull, errMicIdx, sweepCore_scaled, cfg)
% 基于扫频起始点对齐
fs = cfg.fs();

% 理论扫频起始位置
sweep_start_theory = round(cfg.padLeading * fs) + 1;

% 使用扫频信号的自相关来精确定位扫频起始点
mainMicChannel = errMicIdx(1);
recorded_main = recordedFull(:, mainMicChannel);

% 在理论位置附近搜索扫频响应的起始点
search_start = max(1, sweep_start_theory - round(0.05 * fs)); % 前后50ms
search_end = min(length(recorded_main), sweep_start_theory + round(0.1 * fs)); % 前后100ms

% 截取搜索范围内的信号
search_segment = recorded_main(search_start:search_end);

% 使用扫频信号的起始部分进行匹配（修正长度问题）
% 修正：使用嵌套min函数
temp_min = min(length(sweepCore_scaled), length(search_segment));
match_len = min(1000, temp_min); % 匹配前1000个样本或更短
sweep_start_part = sweepCore_scaled(1:match_len);
search_segment_matched = search_segment(1:match_len); % 截取相同长度

% 计算互相关
if length(search_segment_matched) >= 100 && length(sweep_start_part) >= 100  % 确保有足够的样本
    correlation = xcorr(search_segment_matched, sweep_start_part, 'none');
    
    % 归一化相关系数
    norm1 = sqrt(sum(search_segment_matched.^2));
    norm2 = sqrt(sum(sweep_start_part.^2));
    if norm1 > 0 && norm2 > 0
        correlation_normalized = correlation / (norm1 * norm2);
    else
        correlation_normalized = correlation;
    end
    
    lags = -(length(search_segment_matched) - 1):(length(search_segment_matched) - 1);
    
    % 找到最大相关性位置
    [max_corr, max_idx] = max(abs(correlation_normalized));
    lag_at_max = lags(max_idx);
    
    % 计算实际扫频起始位置
    actual_start_idx = sweep_start_theory + lag_at_max;
    
    % 确保在有效范围内
    actual_start_idx = max(search_start, min(search_end, actual_start_idx));
    
    fprintf('    扫频起始点对齐: 理论@%d, 实际@%d (相关性: %.3f)\n', ...
        sweep_start_theory, actual_start_idx, max_corr);
else
    % 如果搜索段太短，使用理论位置
    actual_start_idx = sweep_start_theory;
    fprintf('    扫频起始点对齐: 使用理论位置@%d (搜索段太短)\n', actual_start_idx);
end

% 计算扫频响应段的结束位置
sweep_end_idx = actual_start_idx + length(sweepCore_scaled) - 1;

if sweep_end_idx > size(recordedFull, 1)
    fprintf('    警告：录制长度不足，补零\n');
    recordedFull(end+1:sweep_end_idx, :) = 0;
end

% 截取对齐后的扫频响应段
recorded_aligned = recordedFull(actual_start_idx:sweep_end_idx, errMicIdx);
end