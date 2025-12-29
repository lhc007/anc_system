function mean_coherence = compute_mean_coherence(recorded, sweep, fs)
    % 修复：确保长度匹配
    L = min(length(recorded), length(sweep));
    if L < 1024
        mean_coherence = 0;
        return;
    end
    
    recorded = recorded(1:L);
    sweep = sweep(1:L);
    
    % 去直流（必要）
    recorded = recorded - mean(recorded);
    sweep = sweep - mean(sweep);
    
    % 确保信号能量足够
    if max(abs(recorded)) < 1e-6 || max(abs(sweep)) < 1e-6
        mean_coherence = 0;
        return;
    end
    
    % 归一化信号
    recorded = recorded / (std(recorded) + eps);
    sweep = sweep / (std(sweep) + eps);
    
    % 不要加全局窗！mscohere 内部处理
    try
        % 使用合理参数：段长2048，重叠50%，避免过长段
        [Cxy, f] = mscohere(recorded, sweep, hann(2048), 1024, [], fs);
        
        % 关注 ANC 频段 (20–1000 Hz) - 扩展频段
        freq_mask = (f >= 20) & (f <= 1000);
        valid_vals = Cxy(freq_mask);
        
        if ~isempty(valid_vals) && all(isfinite(valid_vals))
            mean_coherence = mean(valid_vals);
        else
            mean_coherence = 0;
        end
        
    catch
        mean_coherence = 0; % 宁可保守
    end
    
    % 限制范围
    mean_coherence = max(0, min(1, mean_coherence));
end