function [actual_start_idx, recorded_aligned] = align_sweep_start(recordedFull, exciteInfo, errMicIdx)
% ALIGN_SWEEP_START v5.0: 使用 exciteInfo.clickPos 作为唯一事实源
%
% 输入:
%   recordedFull    : 完整录制信号 [N x M]
%   exciteInfo      : 来自 generate_excitation_signal 的 info 结构体
%   errMicIdx       : 误差麦克风通道索引向量
%
% 输出:
%   actual_start_idx: 对齐后的扫频起始样本索引（绝对位置）
%   recorded_aligned: 对齐后的录制段 [L x numErrMics]

    % 初始化输出（防止未赋值错误）
    actual_start_idx = [];
    recorded_aligned = [];

    try
        fs = exciteInfo.fs;
        mainChan = errMicIdx(1);
        rec_main = recordedFull(:, mainChan);

        % === 直接使用生成时的真实值（不再从 cfg 计算！）===
        theory_click_pos = exciteInfo.clickPos;          % ← 关键：来自激励生成
        theory_sweep_start = exciteInfo.sweepStartIdx;   % ← 同样来自激励生成
        sweep_core_len = exciteInfo.coreLength;
        padLeadSamples = exciteInfo.padLeadingSamples;
        padTrailSamples = exciteInfo.padTrailingSamples;

        fprintf('    理论 Click 位置: %d (%.3f秒)\n', theory_click_pos, theory_click_pos/fs);
        fprintf('    理论扫频起始: %d (%.3f秒)\n', theory_sweep_start, theory_sweep_start/fs);

        if theory_click_pos < 1 || theory_click_pos > length(rec_main)
            error('Click 位置超出录制范围');
        end

        % === 搜索 Click ===
        search_win = round(0.5 * fs);  % ±0.5秒窗口
        start_search = max(1, theory_click_pos - search_win);
        end_search = min(length(rec_main), theory_click_pos + search_win);

        if start_search > end_search
            error('Click 搜索窗口无效');
        end

        [~, idx_local] = max(abs(rec_main(start_search:end_search)));
        actual_click_pos = start_search + idx_local - 1;
        system_offset = actual_click_pos - theory_click_pos;
        corrected_sweep_start = theory_sweep_start + system_offset;

        fprintf('    检测到 Click @ %d (%.3f秒), 系统偏移: %+d 样本 (%+.3f秒)\n', ...
            actual_click_pos, actual_click_pos/fs, system_offset, system_offset/fs);

        % === 提取对齐后的扫频段 ===
        end_idx = corrected_sweep_start + sweep_core_len - 1;

        % === 边界处理 ===
        if corrected_sweep_start < 1
            warning('校正后起始位置 < 1，截断到 1');
            corrected_sweep_start = 1;
            end_idx = corrected_sweep_start + sweep_core_len - 1;
        end

        if end_idx > size(recordedFull, 1)
            warning('录制数据不足（需要 %d，实际 %d），末尾补零', end_idx, size(recordedFull,1));
            recorded_part = recordedFull(corrected_sweep_start:end, errMicIdx);
            missing = end_idx - size(recordedFull, 1);
            recorded_aligned = [recorded_part; zeros(missing, size(recorded_part,2))];
        else
            recorded_aligned = recordedFull(corrected_sweep_start:end_idx, errMicIdx);
        end

        actual_start_idx = corrected_sweep_start;

        fprintf('    对齐完成: 起始=%d, 长度=%d, 持续时间=%.3f秒\n\n', ...
            actual_start_idx, size(recorded_aligned,1), size(recorded_aligned,1)/fs);

    catch ME
        % 错误回退：确保输出变量有值
        fprintf('    [align] 错误: %s\n', ME.message);
        if isempty(actual_start_idx)
            actual_start_idx = exciteInfo.sweepStartIdx; % 回退到理论扫频起始
        end
        if isempty(recorded_aligned)
            sweep_core_len = exciteInfo.coreLength;
            end_idx = actual_start_idx + sweep_core_len - 1;
            if end_idx <= size(recordedFull,1)
                recorded_aligned = recordedFull(actual_start_idx:end_idx, errMicIdx);
            else
                recorded_aligned = zeros(sweep_core_len, length(errMicIdx));
                warning('对齐失败，使用零信号');
            end
        end
    end
end