function irResult = deconv_ess(recorded, invFilter, cfg, clickAbsPos, sweepCoreLen)
% 必须以 clickAbsPos 为锚点！
    Lh = cfg.irMaxLen; % e.g., 4096
    h_full = conv(recorded(:), invFilter(:), 'full');
    
    % === 关键：IR 主峰应在 clickAbsPos 附近 ===
    search_start = clickAbsPos + 10;
    search_end   = clickAbsPos + min(2000, length(h_full)); % 最多搜 42ms
    
    if search_end > length(h_full), search_end = length(h_full); end
    
    [~, offset] = max(abs(h_full(search_start:search_end)));
    actual_peak = search_start + offset - 1;
    
    % 截取 IR：保留 peak 前 preDelay 样本（如 128），后补至 Lh
    preDelay = getCfgField(cfg, 'deconvPreDelayKeep', 128); % ← 至少 128！
    ir_start = max(1, actual_peak - preDelay);
    ir_end   = ir_start + Lh - 1;
    
    if ir_end > length(h_full)
        irSeg = [h_full(ir_start:end); zeros(ir_end - length(h_full), 1)];
    else
        irSeg = h_full(ir_start:ir_end);
    end
    
    % SNR 估计：用主峰附近 vs 尾部噪声
    sig_idx = max(1, preDelay - 20) : min(Lh, preDelay + 20);
    noise_idx = (Lh - 100) : Lh;
    snr = 10*log10( mean(irSeg(sig_idx).^2) / (mean(irSeg(noise_idx).^2) + eps) );
    
    irResult = struct(...
        'ir', irSeg, ...
        'snr', snr, ...
        'peakInSegment', preDelay ... % ← 这才是 IR 段内峰值位置（用于 align 和 assess）
    );
end