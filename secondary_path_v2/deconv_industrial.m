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

    % =================== 使用 N + irMaxLen - 1 作为线性卷积长度下限 ===================
    Nfft = 2^nextpow2(N + irMaxLen - 1);

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

    % =================== 修复：安全的IR裁剪 ===================
    % 确定起始位置
    ir_start = max(1, max_idx - preDelayKeep);
    
    % 确保不会越界
    if ir_start > length(h_full)
        % 起始位置超出h_full长度，返回全零IR
        ir = zeros(irMaxLen, 1);
        warning('deconv_industrial: ir_start超出范围，返回全零IR');
        peakPos = 1;  % 设置为默认值
    else
        % 计算结束位置
        ir_end = ir_start + irMaxLen - 1;
        
        if ir_end > length(h_full)
            % 如果结束位置超出，截取可用部分并补零
            available_part = h_full(ir_start:end);
            padding_needed = irMaxLen - length(available_part);
            ir = [available_part; zeros(padding_needed, 1)];
        else
            % 正常情况：直接截取
            ir = h_full(ir_start:ir_end);
        end
        
        % 计算IR内部的峰值位置
        peakPos = max_idx - ir_start + 1;
    end

    % SNR评估
    if peakPos <= length(ir)
        win = max(1, round(0.001*fs));
        sig_start = max(1, peakPos - win);
        sig_end = min(length(ir), peakPos + win);
        noise_end = max(1, sig_start - 1);
        
        if noise_end > 0 && sig_end >= sig_start
            sig_win = ir(sig_start:sig_end);
            noise_win = ir(1:noise_end);
            snr_db = 10*log10(sum(sig_win.^2)/(sum(noise_win.^2)+eps));
        else
            snr_db = -inf;
        end
    else
        % 如果peakPos无效
        snr_db = -inf;
    end

    % 返回结果
    result = struct(...
        'ir', ir, ...
        'peakPos', peakPos, ...
        'snr', snr_db);
end