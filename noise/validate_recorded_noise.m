function validate_recorded_noise()
    filename = 'anc_sim_data_road_noise.wav';
    
    if ~isfile(filename)
        error('❌ 文件不存在: %s', filename);
    end
    
    [x, fs] = audioread(filename);
    [N, C] = size(x);
    
    fprintf('✅ 文件加载成功\n');
    fprintf('   采样率: %d Hz\n', fs);
    fprintf('   时长: %.2f 秒\n', N / fs);
    fprintf('   通道数: %d\n', C);
    
    if C ~= 6
        warning('⚠️ 通道数不是6！当前为 %d', C);
    end
    
    % === 1. 检查是否静音或能量过低 ===
    rms_all = sqrt(mean(x.^2));
    fprintf('\n📊 通道RMS能量:\n');
    for ch = 1:C
        fprintf('   Ch%d: %.4f\n', ch, rms_all(ch));
    end
    
    if any(rms_all < 1e-4)
        warning('⚠️ 某些通道能量极低（可能未连接或静音）');
    end
    
    % === 2. 检查是否所有通道相同（复制）===
    all_same = true;
    for ch = 2:C
        if ~isequal(x(:,1), x(:,ch))
            all_same = false;
            break;
        end
    end
    if all_same
        error('❌ 所有通道完全相同！这是软件复制，不可用！');
    else
        fprintf('\n✅ 通道存在差异（非复制）\n');
    end
    
    % === 3. 计算参考麦 vs 误差麦的互相关延迟 ===
    ref_chs = [1 2 3 4];
    err_chs = [5 6];
    
    fprintf('\n⏱️ 检查参考麦是否领先误差麦（理想：正延迟）:\n');
    delay_ok = false;
    max_lag_samples = round(0.01 * fs); % ±10 ms window
    
    for r = ref_chs
        for e = err_chs
            L = min(2*fs, N);
            % 🔧 修复点：xcorr(Err, Ref) → 得到 Err 相对于 Ref 的延迟
            xc = xcorr(x(1:L, e), x(1:L, r), max_lag_samples, 'coeff');
            [~, idx] = max(abs(xc));
            lag_samples = idx - (length(xc) + 1) / 2;  % 正值 = Err 滞后
            delay_ms = lag_samples / fs * 1000;
            
            fprintf('   Ref%d → Err%d: 延迟 = %.1f ms\n', r, e, delay_ms);
            
            if delay_ms > 0.1
                delay_ok = true;
            end
        end
    end
    
    if delay_ok
        fprintf('✅ 至少一对参考→误差存在正向延迟（物理合理）\n');
    else
        warning('⚠️ 所有参考→误差延迟 ≤ 0！可能噪声源位置错误或扬声器未关闭');
    end
    
    % === 4. 检查是否削波 ===
    max_val = max(abs(x(:)));
    if max_val >= 0.99
        warning('⚠️ 信号接近满幅（max=%.3f），可能存在削波', max_val);
    else
        fprintf('\n✅ 无削波（最大幅度: %.3f）\n', max_val);
    end
    
    % === 5. 可视化（可选）===
    figure('Name', '录制数据验证', 'NumberTitle', 'off');
    subplot(2,1,1);
    plot((0:min(5000,N)-1)/fs, x(1:min(5000,N), :));
    title('前5000样本：各通道波形');
    xlabel('时间 (s)');
    legend('Ref1','Ref2','Ref3','Ref4','Err5','Err6', 'Location', 'best');
    grid on;
    
    subplot(2,1,2);
    corr_mat = corrcoef(x(1:min(10000,N), :)');
    imagesc(corr_mat);
    colorbar;
    title('通道间相关系数矩阵（越白越相似）');
    xlabel('通道'); ylabel('通道');
    set(gca, 'XTick', 1:6, 'YTick', 1:6);
    
    fprintf('\n🎉 验证完成！\n');
    if ~all_same && delay_ok && max_val < 0.99
        fprintf('🟢 数据可用！可用于ANC仿真。\n');
    else
        fprintf('🔴 数据存在问题，请检查录制过程。\n');
    end
end