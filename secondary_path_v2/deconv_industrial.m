function result = deconv_industrial(recorded, sweep, cfg)
    fs = cfg.fs;
    regEps = cfg.deconvRegEps;
    irMaxLen = cfg.irMaxLen;
    preDelayKeep = cfg.deconvPreDelayKeep;

    % 确保列向量
    rec = recorded(:);
    exc = sweep(:);

    % 截断到相同长度
    N = min(length(rec), length(exc));
    rec = rec(1:N);
    exc = exc(1:N);

    fprintf('[DEBUG] deconv输入: rec=%d, exc=%d\n', N, N);

    % =================== 频域反卷积 ===================
    % 使用足够长的 FFT 避免混叠
    Nfft = 2^nextpow2(N + irMaxLen - 1);  % 长度应 >= N + irMaxLen - 1

    % 归一化激励信号
    exc_norm = exc / norm(exc);
    EXC_norm = fft(exc_norm, Nfft);
    REC = fft(rec, Nfft);

    magEXC2 = abs(EXC_norm).^2;
    noise_floor = median(abs(REC)) * 0.1;
    regTerm = regEps * max(magEXC2) + noise_floor^2;

    % 维纳滤波
    Hf = (REC .* conj(EXC_norm)) ./ (magEXC2 + regTerm);

    % IFFT 得到线性卷积结果
    h_full = real(ifft(Hf));

    % 去直流
    h_full = h_full - mean(h_full);

    % =================== 寻找主峰 ===================
    h_abs = abs(h_full);
    search_start = round(0.1 * Nfft);
    search_end = round(0.9 * Nfft);
    [max_val, max_idx] = max(h_abs(search_start:search_end));
    max_idx = max_idx + search_start - 1;

    if max_val < 0.1 * max(h_abs)
        max_idx = round(Nfft/2);
    end

    % =================== SNR计算 ===================
    win_samples = round(0.0005 * fs);
    sig_start = max(1, max_idx - win_samples);
    sig_end = min(Nfft, max_idx + win_samples);

    noise_start = 1;
    noise_end = max(1, sig_start - 1);

    signal_energy = sum(h_full(sig_start:sig_end).^2);
    noise_energy = sum(h_full(noise_start:noise_end).^2);

    if noise_end > noise_start
        noise_energy = noise_energy * (sig_end-sig_start+1) / (noise_end-noise_start+1);
    else
        noise_energy = 1e-6;
    end

    snr_db = 10*log10(signal_energy / noise_energy);

    % =================== IR裁剪 ===================
    ir_start = max(1, max_idx - preDelayKeep);
    ir_end = min(Nfft, ir_start + irMaxLen - 1);
    ir = h_full(ir_start:ir_end);

    if length(ir) < irMaxLen
        ir = [ir; zeros(irMaxLen - length(ir), 1)];
    end

    peak_pos = max_idx - ir_start + 1;
    if peak_pos < 1
        peak_pos = 1;
    elseif peak_pos > irMaxLen
        peak_pos = irMaxLen;
    end

    % =================== 重建验证 ===================
    reconstructed = conv(ir, exc, 'same');
    recon_length = min(length(rec), length(reconstructed));

    mse = mean((rec(1:recon_length) - reconstructed(1:recon_length)).^2);
    recon_snr = 10*log10(sum(rec(1:recon_length).^2)/(mse*recon_length + 1e-12));

    fprintf('[DEBUG] IR提取结果: 长度=%d, SNR=%.1f dB, 重建SNR=%.1f dB\n', ...
        length(ir), snr_db, recon_snr);

    result = struct(...
        'ir', ir, ...
        'peakPos', peak_pos, ...
        'snr', snr_db, ...
        'reconstructed', reconstructed(1:recon_length), ...
        'reconSNR', recon_snr);
end