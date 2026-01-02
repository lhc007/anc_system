function [actual_start_idx, recorded_aligned] = align_sweep_start(recordedFull, sweepSig, cfg, errMicIdx)
% 对齐扫频起始点，并返回对齐后的误差麦克风信号

fs = cfg.fs;
padLeading = cfg.padLeading;
padTrailing = cfg.padTrailing;

% 理论扫频开始位置
sweep_start_theory = round(padLeading * fs) + 1;

% 计算核心扫频长度（总长减去前后静音）
total_len = length(sweepSig);
leading_samples = round(padLeading * fs);
trailing_samples = round(padTrailing * fs);
sweep_core_len = total_len - leading_samples - trailing_samples;

% 使用第一个误差麦克风通道进行对齐
mainChan = errMicIdx(1);
rec_main = recordedFull(:, mainChan);

% 在理论位置附近搜索（±100ms）
search_win = round(0.1 * fs);
start_search = max(1, sweep_start_theory - search_win);
end_search = min(length(rec_main), sweep_start_theory + search_win);
search_seg = rec_main(start_search:end_search);

match_len = min(1000, length(search_seg));
if match_len < 100
    actual_start_idx = sweep_start_theory;
    fprintf('    [align] 搜索段太短，使用理论起始点\n');
else
    ref_part = sweepSig(sweep_start_theory : sweep_start_theory + match_len - 1);
    corr = xcorr(search_seg(1:match_len), ref_part, 'none');
    [~, max_idx] = max(abs(corr));
    lag = (length(ref_part) - 1) - (max_idx - 1);
    actual_start_idx = start_search + lag;
    actual_start_idx = max(1, min(length(rec_main), actual_start_idx));
end

% 截取对齐后的扫频响应段
end_idx = actual_start_idx + sweep_core_len - 1;

if end_idx > size(recordedFull, 1)
    % 补零到所需长度
    padding_rows = end_idx - size(recordedFull, 1);
    recordedFull = [recordedFull; zeros(padding_rows, size(recordedFull, 2))];
end

recorded_aligned = recordedFull(actual_start_idx:end_idx, errMicIdx);

fprintf('    对齐后信号长度: %d (期望: %d)\n', ...
    size(recorded_aligned, 1), sweep_core_len);
end