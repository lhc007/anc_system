function [actual_start_idx, recorded_aligned] = align_sweep_start(recordedFull, sweepSig, cfg, errMicIdx)
% 修复：不修改原始数据的安全对齐函数
% 版本：v2.0 - 增强安全性和健壮性
% 输入：
%   recordedFull - 完整录制的数据 [samples × channels]
%   sweepSig    - 理论扫频信号
%   cfg         - 配置结构体，包含fs, padLeading, padTrailing
%   errMicIdx   - 误差麦克风通道索引
% 输出：
%   actual_start_idx - 实际找到的起始位置
%   recorded_aligned - 对齐后的录制数据

% 参数提取
fs = cfg.fs;
padLeading = cfg.padLeading;
padTrailing = cfg.padTrailing;

% 理论扫频开始位置
sweep_start_theory = round(padLeading * fs) + 1;
fprintf('    理论起始位置: %d (%.3f秒)\n', sweep_start_theory, padLeading);

% 计算核心扫频长度
total_len = length(sweepSig);
leading_samples = round(padLeading * fs);
trailing_samples = round(padTrailing * fs);
sweep_core_len = total_len - leading_samples - trailing_samples;
fprintf('    核心扫频长度: %d 样本 (%.3f秒)\n', sweep_core_len, sweep_core_len/fs);

% 使用第一个误差麦克风通道进行对齐
mainChan = errMicIdx(1);
rec_main = recordedFull(:, mainChan);

% 在理论位置附近搜索（±100ms）
search_win = round(0.1 * fs);
start_search = max(1, sweep_start_theory - search_win);
end_search = min(length(rec_main), sweep_start_theory + search_win);

% 验证搜索范围
if start_search >= end_search
    fprintf('    [align] 警告: 搜索范围无效，使用理论起始点\n');
    actual_start_idx = sweep_start_theory;
else
    search_seg = rec_main(start_search:end_search);
    fprintf('    搜索范围: %d-%d (共%d样本)\n', start_search, end_search, length(search_seg));
    
    % 确定匹配长度（安全边界）
    match_len = min(1000, length(search_seg));
    if match_len < 50
        fprintf('    [align] 搜索段太短，使用理论起始点\n');
        actual_start_idx = sweep_start_theory;
    else
        % 安全截取参考信号（防止越界）
        ref_end = min(sweep_start_theory + match_len - 1, length(sweepSig));
        ref_part = sweepSig(sweep_start_theory : ref_end);
        actual_match_len = length(ref_part);
        
        % 确保搜索信号有足够长度
        if actual_match_len > length(search_seg)
            fprintf('    [align] 警告: 调整匹配长度 %d -> %d\n', ...
                actual_match_len, length(search_seg));
            actual_match_len = length(search_seg);
            % 重新截取参考信号
            ref_end = min(sweep_start_theory + actual_match_len - 1, length(sweepSig));
            ref_part = sweepSig(sweep_start_theory : ref_end);
        end
        
        fprintf('    匹配长度: %d 样本\n', actual_match_len);
        
        % 计算互相关
        corr = xcorr(search_seg(1:actual_match_len), ref_part, 'none');
        
        % 找到最大相关位置
        [~, max_idx] = max(abs(corr));
        lag = (actual_match_len - 1) - (max_idx - 1);
        actual_start_idx = start_search + lag;
        
        % 确保在有效范围内
        actual_start_idx = max(1, min(length(rec_main), actual_start_idx));
        
        % 验证对齐位置是否合理
        offset = actual_start_idx - sweep_start_theory;
        if abs(offset) > search_win
            fprintf('    [align] 警告: 对齐偏移过大 (%d样本, %.3f秒)\n', ...
                offset, offset/fs);
            fprintf('    检查信号质量，可能需要重新录制\n');
        else
            fprintf('    对齐偏移: %d 样本 (%.3f秒)\n', offset, offset/fs);
        end
    end
end

% 安全截取（避免越界）
end_idx = actual_start_idx + sweep_core_len - 1;
actual_end_idx = min(size(recordedFull, 1), end_idx);

% 验证截取范围
if actual_end_idx < actual_start_idx
    error('对齐失败：起始位置 %d 大于结束位置 %d', actual_start_idx, actual_end_idx);
end

if actual_end_idx > size(recordedFull, 1)
    fprintf('    [align] 警告: 结束位置 %d 超出数据长度 %d\n', ...
        actual_end_idx, size(recordedFull, 1));
end

% 截取对齐段
extracted_len = actual_end_idx - actual_start_idx + 1;
recorded_aligned = recordedFull(actual_start_idx:actual_end_idx, errMicIdx);

% 如果截取长度不足，补零
if extracted_len < sweep_core_len
    padding_needed = sweep_core_len - extracted_len;
    recorded_aligned = [recorded_aligned; zeros(padding_needed, size(recorded_aligned, 2))];
    fprintf('    [align] 警告: 截取长度不足，补零 %d 样本 (%.3f秒)\n', ...
        padding_needed, padding_needed/fs);
elseif extracted_len > sweep_core_len
    fprintf('    [align] 注意: 截取长度 %d 超过期望 %d\n', ...
        extracted_len, sweep_core_len);
end

% 最终验证
fprintf('    对齐完成: 起始位置 %d, 长度 %d (期望: %d)\n', ...
    actual_start_idx, size(recorded_aligned, 1), sweep_core_len);
fprintf('    持续时间: %.3f秒 (期望: %.3f秒)\n\n', ...
    size(recorded_aligned, 1)/fs, sweep_core_len/fs);

% 可选：返回实际截取长度信息（如果需要）
% if nargout > 2
%     varargout{1} = extracted_len;
% end
end