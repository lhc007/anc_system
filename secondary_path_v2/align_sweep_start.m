function [actual_start_idx, recorded_aligned, clickAbsPos] = align_sweep_start(recordedFull, exciteInfo, errMicIdx)
    % ALIGN_SWEEP_START v6.0: 使用 Sweep 能量起始对齐（无 Click）

    actual_start_idx = [];
    recorded_aligned = [];
    clickAbsPos = [];

    try
        fs = exciteInfo.fs;
        mainChan = errMicIdx(1);
        rec_main = recordedFull(:, mainChan);

        % === 理论扫频起始位置（来自激励生成）===
        theory_sweep_start = exciteInfo.sweepStartIdx;
        sweep_core_len = exciteInfo.coreLength;

        fprintf(' 理论扫频起始: %d (%.3f秒)\n', theory_sweep_start, theory_sweep_start/fs);

        % === 搜索窗口：围绕理论起始点 ±0.5秒 ===
        search_win = round(0.5 * fs);
        start_search = max(1, theory_sweep_start - search_win);
        end_search = min(length(rec_main), theory_sweep_start + search_win);

        if start_search >= end_search
            error('搜索窗口无效');
        end

        window_sig = rec_main(start_search:end_search);

        % === 方法：使用能量导数找上升沿 ===
        energy = abs(window_sig).^2;
        smooth_energy = movmean(energy, max(1, round(0.005 * fs))); % 5ms 平滑
        d_energy = diff([0; smooth_energy]); % 一阶导数

        % 找最大正跳变（能量最快上升点）
        [~, rel_peak] = max(d_energy);

        % 验证：确保该点确实有显著能量
        if smooth_energy(rel_peak) < 0.01 * max(smooth_energy)
            warning('未检测到显著能量上升，回退到理论位置');
            actual_sweep_start = theory_sweep_start;
        else
            actual_sweep_start = start_search + rel_peak - 1;
        end

        system_offset = actual_sweep_start - theory_sweep_start;
        fprintf(' 检测到 Sweep 起始 @ %d (%.3f秒), 系统偏移: %+d 样本 (%+.3f秒)\n', ...
            actual_sweep_start, actual_sweep_start/fs, system_offset, system_offset/fs);

        % === 提取对齐后的扫频段 ===
        end_idx = actual_sweep_start + sweep_core_len - 1;

        if actual_sweep_start < 1
            warning('校正后起始位置 < 1，截断到 1');
            actual_sweep_start = 1;
            end_idx = actual_sweep_start + sweep_core_len - 1;
        end

        if end_idx > size(recordedFull, 1)
            warning('录制数据不足（需要 %d，实际 %d），末尾补零', end_idx, size(recordedFull,1));
            recorded_part = recordedFull(actual_sweep_start:end, errMicIdx);
            missing = end_idx - size(recordedFull, 1);
            recorded_aligned = [recorded_part; zeros(missing, size(recorded_part,2))];
        else
            recorded_aligned = recordedFull(actual_sweep_start:end_idx, errMicIdx);
        end

        actual_start_idx = actual_sweep_start;
        clickAbsPos = actual_sweep_start; % ⬅️ 关键：传递真实 sweep 起始位置给 deconv_ess

        fprintf(' 对齐完成: 起始=%d, 长度=%d, 持续时间=%.3f秒\n\n', ...
            actual_start_idx, size(recorded_aligned,1), size(recorded_aligned,1)/fs);

    catch ME
        fprintf(' [align] 错误: %s\n', ME.message);
        % 回退到理论值
        actual_start_idx = exciteInfo.sweepStartIdx;
        clickAbsPos = actual_start_idx;
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