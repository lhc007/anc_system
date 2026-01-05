function irResult = deconv_ess(recorded, invFilter, cfg)
%DECONV_ESS 使用ESS反卷积从录制信号中提取冲激响应
%
% 输入:
%   recorded    : 对齐后的录制信号（单通道，列向量）
%   invFilter   : 反卷积滤波器（= flipud(weightedSweepCore)）
%   cfg         : 配置结构体，需包含 irMaxLen
%
% 输出:
%   irResult.ir       : 提取的冲激响应 [Lh x 1]
%   irResult.snr      : 估计的SNR (dB)
%   irResult.peakPos  : 主峰在 ir 中的位置（样本索引，1-based）

    % --- 参数 ---
    Lh = cfg.irMaxLen;
    if Lh <= 0
        error('irMaxLen 必须为正整数');
    end

    % 确保输入为列向量
    recorded = recorded(:);
    invFilter = invFilter(:);

    % === 标准 ESS 反卷积：使用 'same' 模式 ===
    % 理论依据：ESS 的自相关主瓣出现在末端，'same' 保留中心区域（含完整 IR）
    h_full = conv(recorded, invFilter, 'same');  % 长度 = length(recorded)

    % === 关键：IR 能量集中在 h_full 的末尾附近 ===
    % Farina (2000) 指出：IR 出现在反卷积结果的最后 ~T 秒内
    N = length(h_full);
    if N < Lh
        % 数据太短，前面补零
        ir = [zeros(Lh - N, 1); h_full];
    else
        % 截取最后 Lh 个样本（确保包含主峰）
        ir = h_full(end - Lh + 1 : end);
    end

    % 再次确保长度正确（防御性编程）
    if length(ir) > Lh
        ir = ir(1:Lh);
    elseif length(ir) < Lh
        ir = [zeros(Lh - length(ir), 1); ir];
    end

    % === 找主峰位置（在 ir 内部）===
    [~, peakIdx] = max(abs(ir));
    
    % === SNR 估计：主峰前 vs 主峰后噪声 ===
    snr_db = estimate_snr_from_tail(ir, peakIdx);

    % === 输出 ===
    irResult.ir = ir;
    irResult.snr = snr_db;
    irResult.peakPos = peakIdx;  % 1-based index in ir
end

% -------------------------------------------------------------------------
function snr_db = estimate_snr_from_tail(ir, peakPos)
% 估计 SNR：主峰能量 vs 尾部噪声能量
% 假设主峰后 50% 区域为噪声（避免预回声干扰）

    L = length(ir);
    if peakPos >= L
        snr_db = -Inf;
        return;
    end

    % 主峰能量（取峰值附近小窗口，更鲁棒）
    win = min(5, floor((L - peakPos)/2));
    signal_energy = mean(abs(ir(peakPos : min(L, peakPos + win))).^2);

    % 噪声区域：主峰之后的后半段（避开可能的反射）
    noise_start = min(L, peakPos + 10);  % 跳过主峰尾迹
    noise_end = L;
    if noise_end <= noise_start
        noise_energy = eps;
    else
        noise_samples = ir(noise_start:noise_end);
        noise_energy = mean(abs(noise_samples).^2);
    end

    if noise_energy <= 0
        noise_energy = eps;
    end

    snr_db = 10 * log10(signal_energy / noise_energy);
end