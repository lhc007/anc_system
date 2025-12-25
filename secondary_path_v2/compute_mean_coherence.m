function mean_coherence = compute_mean_coherence(recorded, sweep, fs)
% 计算记录信号与激励信号之间的平均相干性
% recorded: 记录信号
% sweep: 激励信号
% fs: 采样率
% mean_coherence: 关键频段的平均相干性

    % 确保长度一致
    L = min(length(recorded), length(sweep));
    recorded = recorded(1:L);
    sweep = sweep(1:L);
    
    % 检查信号是否非零
    if all(recorded == 0) || all(sweep == 0)
        mean_coherence = 0;
        return;
    end
    
    % 计算相干性
    window_size = min(1024, floor(L/8)); % 自适应窗口大小
    noverlap = floor(window_size/2);
    nfft = max(256, 2^nextpow2(window_size));
    
    try
        [Cxy, f] = mscohere(recorded, sweep, hann(window_size), noverlap, nfft, fs);
        
        % 计算关键频段（ANC关注频段20-1000Hz）的平均相干性
        freq_range = f >= 20 & f <= 1000;
        if any(freq_range)
            mean_coherence = mean(Cxy(freq_range));
        else
            mean_coherence = 0;
        end
    catch
        % 如果mscohere失败，返回保守估计
        mean_coherence = 0;
    end
end