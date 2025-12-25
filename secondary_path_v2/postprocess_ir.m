function irOut = postprocess_ir(irIn, cfg)
% IR后处理函数
% 专门为ANC次级路径测量设计的后处理
% 包括去直流、加窗和可选滤波，但不进行归一化（保留物理增益）
% 输入:
%   irIn: [numSamples, numMics] 输入的脉冲响应
%   cfg: 配置结构体
% 输出:
%   irOut: [numSamples, numMics] 处理后的脉冲响应

[numSamples, numMics] = size(irIn);

% 检查输入维度
if numMics == 1 && size(irIn, 1) > size(irIn, 2)
    % 单通道情况，确保列向量
    irIn = irIn(:);
    [numSamples, numMics] = size(irIn);
end

irOut = zeros(size(irIn));

for m = 1:numMics
    % 提取当前通道的IR
    ir = irIn(:, m);
    
    % ============== 1. 去直流偏置 ==============
    % 使用更稳健的去直流方法
    % 先找到主要信号区域，避免用噪声区域估计直流
    if length(ir) > 100
        % 找到信号峰值
        [~, peakIdx] = max(abs(ir));
        
        % 选择无信号区域估计直流（峰值前50个样本）
        noiseStart = max(1, peakIdx - 100);
        noiseEnd = max(1, peakIdx - 50);
        
        if noiseEnd > noiseStart
            dcEstimate = mean(ir(noiseStart:noiseEnd));
        else
            dcEstimate = mean(ir(1:min(50, length(ir))));
        end
    else
        % 短IR，使用全部样本
        dcEstimate = mean(ir);
    end
    
    % 减去直流偏置
    ir = ir - dcEstimate;
    
    % ============== 2. 加窗处理（可选） ==============
    % 减少边界效应和截断噪声
    if isfield(cfg, 'applyWindow') && cfg.applyWindow
        % 创建自适应窗函数
        % 找到IR的有效长度（基于能量衰减）
        ir_abs = abs(ir);
        energy_cum = cumsum(ir_abs.^2);
        total_energy = energy_cum(end);
        
        if total_energy > 0
            % 找到能量达到95%的位置
            energy_threshold = 0.95 * total_energy;
            valid_end = find(energy_cum >= energy_threshold, 1, 'first');
            
            if isempty(valid_end) || valid_end < numSamples * 0.1
                valid_end = round(numSamples * 0.8); % 默认80%
            end
            
            % 创建Tukey窗，在有效长度后逐渐衰减
            win = ones(numSamples, 1);
            fade_start = max(valid_end, round(numSamples * 0.5)); % 从50%或有效结束点开始
            fade_len = numSamples - fade_start + 1;
            
            if fade_len > 10
                fade_win = tukeywin(fade_len * 2, 0.5);
                win(fade_start:end) = fade_win(fade_len+1:end);
            end
            
            ir = ir .* win;
        else
            % 零能量，不进行加窗
        end
    end
    
    % ============== 3. 频带限制滤波（可选） ==============
    % 限制到ANC系统关心的频带，减少带外噪声
    if isfield(cfg, 'applyPostFilter') && cfg.applyPostFilter && cfg.fs > 1000
        % ANC通常关注20-1000Hz频段
        f_low = 20;  % 下限频率
        f_high = min(1000, 0.45 * cfg.fs);  % 上限频率，考虑抗混叠
        
        % 使用线性相位FIR滤波器，保持相位特性
        filt_order = 100;  % 适当阶数
        b = fir1(filt_order, [f_low, f_high] / (cfg.fs/2), 'bandpass');
        
        % 零相位滤波
        ir = filtfilt(b, 1, ir);
    end
    
    % ============== 4. 轻微噪声抑制（可选） ==============
    % 抑制尾部噪声，减少预回声
    if isfield(cfg, 'applyNoiseGate') && cfg.applyNoiseGate
        % 找到信号峰值
        [peakVal, peakIdx] = max(abs(ir));
        
        % 估计噪声水平（使用前10%的数据）
        noise_start = 1;
        noise_end = max(1, min(round(numSamples * 0.1), peakIdx - 10));
        
        if noise_end > noise_start
            noise_segment = abs(ir(noise_start:noise_end));
            noise_level = median(noise_segment);
            
            % 设置门限（-40dB或噪声水平的100倍）
            gate_threshold = max(noise_level * 100, peakVal * 0.01);
            
            % 应用软门限
            for i = 1:length(ir)
                if abs(ir(i)) < gate_threshold
                    % 软衰减
                    attenuation = (abs(ir(i)) / gate_threshold)^2;
                    ir(i) = ir(i) * attenuation;
                end
            end
        end
    end
    
    % 存储处理后的IR
    irOut(:, m) = ir;
end
end