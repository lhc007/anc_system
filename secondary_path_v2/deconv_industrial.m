function result = deconv_industrial(recorded, sweep, cfg)
% 反卷积：从录制信号中提取冲激响应（IR）
% 修复点：正确计算FFT长度，避免时域混叠

    fs = cfg.fs;
    regEps = cfg.deconvRegEps;
    irMaxLen = cfg.irMaxLen;
    preDelayKeep = cfg.deconvPreDelayKeep;

    rec = recorded(:);
    exc = sweep(:);
    N = min(length(rec), length(exc));
    rec = rec(1:N);
    exc = exc(1:N);

    % =================== 修正：使用 N + irMaxLen - 1 作为线性卷积长度下限 ===================
    Nfft = 2^nextpow2(N + irMaxLen - 1);  % ←←← 关键修复！

    % 归一化激励
    exc_norm = exc / norm(exc);
    EXC = fft(exc_norm, Nfft);
    REC = fft(rec, Nfft);

    % 正则化维纳反卷积
    magEXC2 = abs(EXC).^2;
    noise_floor = median(abs(REC(1:round(N/10)))); % 前10%估计噪声
    regTerm = regEps * max(magEXC2) + (noise_floor * 0.1)^2;
    Hf = (REC .* conj(EXC)) ./ (magEXC2 + regTerm);
    h_full = real(ifft(Hf));

    % 寻找主峰（限制在物理合理范围）
    search_start = max(1, round(cfg.minPhysDelaySamples * 0.5));
    search_end = min(length(h_full), round(cfg.maxPhysDelaySamples * 1.5));
    if search_end > search_start
        [~, max_idx] = max(abs(h_full(search_start:search_end)));
        max_idx = max_idx + search_start - 1;
    else
        max_idx = round(length(h_full)/2);
    end

    % 裁剪IR（保留preDelayKeep个前置样本）
    ir_start = max(1, max_idx - preDelayKeep);
    ir_end = ir_start + irMaxLen - 1;
    if ir_end > length(h_full)
        ir = [h_full(ir_start:end); zeros(irEnd - length(h_full), 1)];
    else
        ir = h_full(ir_start:ir_end);
    end

    peakPos = max_idx - ir_start + 1; % IR内部峰值位置（从1开始）

    % SNR评估
    win = max(1, round(0.001*fs));
    sig_win = ir(max(1,peakPos-win):min(end,peakPos+win));
    noise_win = ir(1:max(1,peakPos-win-1));
    snr_db = 10*log10(sum(sig_win.^2)/(sum(noise_win.^2)+eps));

    % 返回结果（包含peakPos供后续使用）
    result = struct(...
        'ir', ir, ...
        'peakPos', peakPos, ...
        'snr', snr_db);
end