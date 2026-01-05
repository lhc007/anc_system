function coh = compute_reconstruction_coherence(ir_extracted, weightedSweep, recorded_segment)
    % 确保列向量
    ir_extracted = ir_extracted(:);
    weightedSweep = weightedSweep(:);
    recorded_segment = recorded_segment(:);

    % 重建信号（'same' 模式 —— 正确！）
    reconstructed = conv(ir_extracted, weightedSweep, 'same');

    % 强制长度对齐（防御性编程）
    L = length(recorded_segment);
    if length(reconstructed) > L
        reconstructed = reconstructed(1:L);
    elseif length(reconstructed) < L
        reconstructed = [reconstructed; zeros(L - length(reconstructed), 1)];  % ← 修复此处
    end

    % 自适应窗口设置（防止短信号崩溃）
    win_len = min(1024, floor(L / 2));
    if win_len < 64
        win_len = min(L, 64);
    end
    window = hamming(win_len);
    noverlap = floor(win_len / 2);
    nfft = max(512, 2^nextpow2(win_len));

    % 计算功率谱和互谱
    Pxy = cpsd(recorded_segment, reconstructed, window, noverlap, nfft);
    Pxx = pwelch(recorded_segment, window, noverlap, nfft);
    Pyy = pwelch(reconstructed, window, noverlap, nfft);

    % 计算相干性（加 eps 防除零）
    coherence = abs(Pxy).^2 ./ (Pxx .* Pyy + eps);

    % 排除不可靠频段（DC 和高频尾部）
    valid_start = 2;               % 跳过 DC (index=1)
    valid_end = length(coherence) - 5;  % 跳过混叠区
    if valid_end <= valid_start
        coh = mean(coherence);  % 信号太短，无法裁剪
    else
        valid_coherence = coherence(valid_start:valid_end);
        % 剔除 NaN/Inf
        valid_coherence = valid_coherence(~isnan(valid_coherence) & ~isinf(valid_coherence));
        if isempty(valid_coherence)
            coh = 0;
        else
            coh = mean(valid_coherence);
        end
    end
end