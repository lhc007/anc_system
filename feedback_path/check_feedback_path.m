%% === 检查反馈路径可用性 ===
function isValid = check_feedback_path(F, delayEst, fs, numRef, numSpk)
    isValid = true;
    fprintf('\n🔍 反馈路径可用性检查:\n');
    
    % 1. 检查是否全零
    totalEnergy = sum(F(:).^2);
    if totalEnergy < 1e-12
        warning('❌ 所有反馈路径能量过低（可能全零）！');
        isValid = false;
        return;
    end
    
    % 2. 逐通道检查
    minDelaySamples = round(0.001 * fs);   % 1ms 最小物理延迟
    maxDelaySamples = round(0.15 * fs);    % 150ms 最大合理延迟（管道+电路）
    
    for s = 1:numSpk
        for r = 1:numRef
            ir = F(:, r, s);
            energy = sum(ir.^2);
            
            % (a) 能量过低？
            if energy < 1e-10
                warning('⚠️  Spk%d→Ref%d: 能量过低（%.2e）', s, r, energy);
                isValid = false;
                continue;
            end
            
            % (b) 主峰是否明显？
            [peakVal, peakIdx] = max(abs(ir));
            tailEnergy = sum(ir(peakIdx+1:end).^2);
            peakRatio = peakVal^2 / (energy + 1e-15);
            if peakRatio < 0.05
                warning('⚠️  Spk%d→Ref%d: 主峰不明显（峰值占比 %.1f%%）', ...
                    s, r, peakRatio*100);
                isValid = false;
            end
            
            % (c) 延迟是否合理？
            d = delayEst(r, s);
            if d < minDelaySamples || d > maxDelaySamples
                warning('⚠️  Spk%d→Ref%d: 延迟异常 (%d 样本 ≈ %.1f ms)', ...
                    s, r, d, 1000*d/fs);
                isValid = false;
            end
        end
    end
    
    % 3. 同一扬声器的延迟一致性（可选）
    for s = 1:numSpk
        delays = delayEst(:, s);
        iqr_val = iqr(delays);
        if iqr_val > 20 % 样本（@48kHz ≈ 0.4ms）
            warning('⚠️  Spk%d: 到各参考麦的延迟差异过大（IQR = %d 样本）', s, iqr_val);
        end
    end
    
    if isValid
        fprintf('✅ 所有反馈路径通过基本可用性检查。\n');
    else
        fprintf('❗ 部分通道存在异常，请检查硬件连接或增益设置。\n');
    end
end