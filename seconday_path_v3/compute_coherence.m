function coh = compute_coherence(y_rec, sweep_core, ir_est, ir_offsets, cfg, varargin)
% COMPUTE_COHERENCE: 基于预测信号对齐的次级路径相干性评估
% 输入:
%   y_rec       : [N × M] 录制信号（多通道）
%   sweep_core  : [L_sweep × 1] 核心激励信号（不含静音）
%   ir_est      : [L_ir × M] 估计的冲激响应
%   ir_offsets  : [1 × M] IR 在原始录制信号中的起始索引（全局坐标）
%   cfg         : 配置结构体（需含 .fs）
%   'Verbose'   : (可选) 是否打印详细对齐与相干性计算过程
%
% 输出:
%   coh         : [1 × M] 各通道的相干性指标 ∈ [0, 1]

    p = inputParser;
    addParameter(p, 'Verbose', false);
    parse(p, varargin{:});
    verbose = p.Results.Verbose;

    fs = cfg.fs;
    numMics = size(y_rec, 2);
    coh = zeros(1, numMics);

    % 分析窗口：200 ms（用于提取局部响应段）
    win_half_samples = round(0.2 * fs); 
    L_sweep = length(sweep_core);

    for ch = 1:numMics
        y_measured = y_rec(:, ch);          % 实测信号
        h_est = ir_est(:, ch);              % 估计IR
        peak_idx = ir_offsets(ch);          % 主峰在 y_measured 中的位置（全局索引）

        % === 1. 生成预测信号 y_pred = sweep_core * h_est ===
        y_pred = conv(sweep_core(:), h_est(:));
        L_pred = length(y_pred);

        % 找预测信号的主峰（更稳健地确定其"中心"）
        [~, pred_peak_idx] = max(abs(y_pred));  % 预测响应的峰值位置

        if verbose
            fprintf('通道 %d | 实测信号长度=%d，预测信号长度=%d\n', ch, length(y_measured), L_pred);
            fprintf('  实测主峰位置：%d（全局索引），预测主峰位置：%d\n', peak_idx, pred_peak_idx);
        end

        % === 2. 提取以主峰为中心的局部窗口 ===
        % 实测信号窗口
        start_meas = max(1, peak_idx - win_half_samples);
        end_meas   = min(length(y_measured), peak_idx + win_half_samples);
        L_meas_win = end_meas - start_meas + 1;

        % 预测信号窗口（以 pred_peak_idx 为中心）
        start_pred = max(1, pred_peak_idx - win_half_samples);
        end_pred   = min(L_pred, pred_peak_idx + win_half_samples);
        L_pred_win = end_pred - start_pred + 1;

        % 取公共长度
        L_common = min(L_meas_win, L_pred_win);

        if verbose
            fprintf('  窗口范围：实测=[%d:%d]（%d点），预测=[%d:%d]（%d点）→ 公共长度=%d\n', ...
                start_meas, end_meas, L_meas_win, start_pred, end_pred, L_pred_win, L_common);
        end

        if L_common > 32
            % 提取等长片段
            y_act = y_measured(start_meas : start_meas + L_common - 1);
            y_recn = y_pred(start_pred : start_pred + L_common - 1);

            % === 3. 互相关精细对齐（±200ms）===
            max_lag_samples = min(round(0.2 * fs), floor(L_common / 2));
            try
                [xc, lags] = xcorr(y_act, y_recn, max_lag_samples, 'none');
                [~, idx_max] = max(abs(xc));
                best_lag = lags(idx_max);  % >0: y_recn 滞后 y_act
            catch
                best_lag = 0;
                if verbose, fprintf('  ⚠️ 互相关计算失败，使用零延迟对齐\n'); end
            end

            if verbose
                fprintf('  最佳对齐延迟：%d 样点（%.3f 秒）\n', best_lag, best_lag / fs);
            end

            % === 4. 对齐后计算归一化相关系数（corr^2）===
            if best_lag > 0
                % y_recn 滞后 → 截断 y_act 前部 或 y_recn 后部
                y_a_aligned = y_act(best_lag + 1 : end);
                y_r_aligned = y_recn(1 : end - best_lag);
            elseif best_lag < 0
                % y_recn 超前
                y_a_aligned = y_act(1 : end + best_lag);
                y_r_aligned = y_recn(-best_lag + 1 : end);
            else
                y_a_aligned = y_act;
                y_r_aligned = y_recn;
            end

            L_aligned = length(y_a_aligned);
            if L_aligned > 16
                % 符号校正（防止反相导致相干性低）
                if sum(y_a_aligned .* y_r_aligned) < 0
                    y_r_aligned = -y_r_aligned;
                end

                % 计算平方相关系数（范围 [0,1]）
                numerator = (y_a_aligned' * y_r_aligned)^2;
                denominator = (sum(y_a_aligned.^2) * sum(y_r_aligned.^2)) + eps;
                coh_val = numerator / denominator;
                coh(ch) = max(0, min(1, coh_val));

                if verbose
                    fprintf('  时域相干性（平方相关系数）：%.6f\n', coh(ch));
                end
            else
                % 对齐后太短 → 回退到频域方法
                if verbose, fprintf('  对齐后信号过短 → 使用频域相干性（mscohere）\n'); end
                coh(ch) = fallback_mscohere(y_act, y_recn, cfg, verbose);
            end
        else
            % 窗口太小 → 直接使用更长片段进行频域相干性估计
            if verbose, fprintf('  局部窗口过短 → 使用频域相干性（mscohere）\n'); end
            coh(ch) = fallback_mscohere(y_measured, y_pred, cfg, verbose);
        end
    end
end

%% -------------------------------------------------------------------------
function c = fallback_mscohere(x, y, cfg, verbose)
% FALLBACK_MSCOHERE: 使用 MATLAB mscohere 计算频域平均相干性
% 当时域对齐不可靠或信号过短时调用

    try
        fs = cfg.fs;
        nfft = 1024;
        win_len = 512;
        noverlap = 256;
        win = hamming(win_len);

        % 确保信号足够长
        if length(x) < win_len || length(y) < win_len
            if verbose, fprintf('  ❌ 信号长度不足，无法计算 mscohere → 相干性设为 0\n'); end
            c = 0;
            return;
        end

        % 计算相干性谱
        [Cxy, F] = mscohere(x, y, win, noverlap, nfft, fs);

        % 选择频带
        if isfield(cfg, 'cohFreqBand') && ~isempty(cfg.cohFreqBand)
            f_band = cfg.cohFreqBand;
            valid_freq = (F >= f_band(1)) & (F <= f_band(2));
            if ~any(valid_freq)
                if verbose, fprintf('  ⚠️ 指定频段 [%g, %g] Hz 超出有效范围 → 改用全频段\n', f_band(1), f_band(2)); end
                valid_freq = (F > 0) & (F <= fs/2);
            end
        else
            valid_freq = (F > 0) & (F <= fs/2);
        end

        c = mean(Cxy(valid_freq));
        c = max(0, min(1, c));  % 安全钳位

        if verbose
            fprintf('  频域相干性（mscohere）：%.6f（基于 %d 个频率点平均）\n', c, sum(valid_freq));
        end

    catch ME
        if verbose
            fprintf('  ❌ mscohere 计算失败：%s → 相干性设为 0\n', ME.message);
        end
        c = 0;
    end
end