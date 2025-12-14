function analyze_speaker_lowfreq()
% analyze_speaker_lowfreq
% 自动分析扬声器低频输出能力（基于 secondary_path.mat）
% 判断是否存在高通滤波（HPF）、低频衰减，并给出 ANC 配置建议

    % --- 加载数据 ---
    if ~exist('secondary_path.mat', 'file')
        error('❌ 文件 secondary_path.mat 不存在，请先运行测量！');
    end
    data = load('secondary_path.mat');
    sec = data.secondary;
    
    fs = sec.fs;
    ir = squeeze(sec.impulseResponses(:,1,1));  % 默认 Spk1, Mic1
    
    if isempty(ir)
        error('❌ 冲激响应为空！');
    end
    
    fprintf('✅ 成功加载 IR，长度 = %d 样本，fs = %d Hz\n', length(ir), fs);
    
    % --- 计算频响 ---
    Nfft = 2^nextpow2(length(ir)*4);  % 零填充提升频率分辨率
    H = fft(ir, Nfft);
    freq = (0:Nfft-1)' * fs / Nfft;
    
    % 关注 20–500 Hz
    idx = freq >= 20 & freq <= 500;
    mag = abs(H(idx));
    mag_db = 20*log10(mag + eps);
    
    % 平滑（1/12 倍频程近似）
    smooth_win = round(0.02 * length(mag_db)); % ~2% 窗长
    if smooth_win < 3, smooth_win = 3; end
    mag_db_smooth = movmean(mag_db, smooth_win);
    
    % --- 关键频率点幅度 ---
    f_ref = 100;  % 参考频率
    [~, i100] = min(abs(freq(idx) - f_ref));
    ref_level = mag_db_smooth(i100);
    
    % 提取关键频点（20, 30, 40, 50, 60, 80 Hz）
    test_freqs = [20, 30, 40, 50, 60, 80];
    levels = zeros(size(test_freqs));
    for k = 1:length(test_freqs)
        [~, ik] = min(abs(freq(idx) - test_freqs(k)));
        levels(k) = mag_db_smooth(ik);
    end
    attenuation = ref_level - levels;  % 相对于 100 Hz 的衰减（dB）
    
    % --- 判断低频能力 ---
    fprintf('\n🔍 低频输出能力分析（相对于 100 Hz）:\n');
    for k = 1:length(test_freqs)
        fprintf('  %3d Hz: %.1f dBFS | 衰减 %.1f dB\n', ...
            test_freqs(k), levels(k), attenuation(k));
    end
    
    % 判断是否存在 HPF
    hpf_suspected = false;
    cutoff_est = NaN;
    
    % 如果 20–40 Hz 衰减 >25 dB，且 60–80 Hz 衰减 <15 dB → 很可能有 HPF
    if attenuation(1) > 25 && attenuation(1) > attenuation(end) + 10
        hpf_suspected = true;
        % 估算截止频率（衰减 = 3 dB 处）
        target_atten = 3;
        if any(attenuation <= target_atten)
            i_cutoff = find(attenuation <= target_atten, 1, 'first');
            cutoff_est = test_freqs(i_cutoff);
        else
            % 找最小衰减点
            [~, i_min] = min(attenuation);
            cutoff_est = test_freqs(i_min);
        end
    end
    
    % --- 判断可用频段 ---
    usable_freq_min = 100; % 默认
    if attenuation(end) <= 12  % 80 Hz 衰减 ≤12 dB → 可用
        usable_freq_min = 80;
    elseif attenuation(5) <= 15  % 60 Hz 衰减 ≤15 dB
        usable_freq_min = 60;
    elseif attenuation(4) <= 18  % 50 Hz
        usable_freq_min = 50;
    else
        usable_freq_min = 100;
    end
    
    % --- 输出结论 ---
    fprintf('\n📋 分析结论:\n');
    if hpf_suspected
        fprintf('  ⚠️  高度怀疑存在高通滤波器（HPF）\n');
        if ~isnan(cutoff_est)
            fprintf('      估算截止频率 ≈ %d Hz\n', cutoff_est);
        end
    else
        fprintf('  ✅ 未检测到明显 HPF，但低频自然衰减严重\n');
    end
    
    fprintf('  📌 建议 ANC 工作频段下限: %d Hz\n', usable_freq_min);
    fprintf('  💡 推荐配置:\n');
    fprintf('      cfg.lowFreqCutHz = %d;\n', max(usable_freq_min - 10, 20));
    fprintf('      cfg.sweepDuration = 16;  %% 提升低频分辨率\n');
    
    % --- 绘图 ---
    figure('Name', '扬声器低频响应分析', 'NumberTitle', 'off');
    semilogx(freq(idx), mag_db_smooth, 'b-', 'LineWidth', 1.5); hold on;
    yline(ref_level, 'k--', '100 Hz 参考');
    xline(usable_freq_min, 'r--', sprintf('建议下限 (%d Hz)', usable_freq_min));
    if hpf_suspected && ~isnan(cutoff_est)
        xline(cutoff_est, 'm--', sprintf('HPF 估算 (%d Hz)', cutoff_est));
    end
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dBFS)');
    title('扬声器电-声频率响应（20–500 Hz）');
    grid on; axis tight;
    xlim([20 500]);
    
    % 标注关键点
    for k = 1:length(test_freqs)
        text(test_freqs(k), levels(k)+1, sprintf('%.0f Hz\n%.1f dB', test_freqs(k), levels(k)), ...
            'HorizontalAlignment','center', 'FontSize',8);
    end
    
    fprintf('\n📊 图形已生成：查看 "扬声器低频响应分析" 窗口\n');
end