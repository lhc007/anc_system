%% 延迟估计验证
function [validatedDelay, result] = validate_delay_estimate(...
    delayEstimate, recordedFull, sweep, cfg)
% % 延迟估计验证
% 检查延迟的物理合理性和信号一致性

fs = cfg.fs;
result = struct();
result.pass = false;
result.confidence = 0;
result.issues = {};

% ============== 检查1: 物理延迟范围 ==============
minPhysical = cfg.minPhysDelaySamples;
maxPhysical = cfg.maxPhysDelaySamples;

if delayEstimate < minPhysical
    result.issues{end+1} = sprintf('延迟(%d)小于最小物理延迟(%d)', ...
        delayEstimate, minPhysical);
end

if delayEstimate > maxPhysical
    result.issues{end+1} = sprintf('延迟(%d)大于最大物理延迟(%d)', ...
        delayEstimate, maxPhysical);
end

% ============== 检查2: 基于声速的物理验证 ==============
% 假设声速343 m/s，计算对应的物理距离
time_sec = delayEstimate / fs;
distance_m = time_sec * 343;

% 合理的物理距离范围（根据应用场景调整）
minDistance = 0.05;  % 5厘米
maxDistance = 5.0;   % 5米

if distance_m < minDistance
    result.issues{end+1} = sprintf('距离(%.3f米)小于最小物理距离(%.3f米)', ...
        distance_m, minDistance);
end

if distance_m > maxDistance
    result.issues{end+1} = sprintf('距离(%.3f米)大于最大物理距离(%.3f米)', ...
        distance_m, maxDistance);
end

% ============== 检查3: 信号验证 ==============
% 提取延迟后的信号段，计算与激励信号的相干性
if delayEstimate + length(sweep) <= length(recordedFull)
    extracted = recordedFull(delayEstimate:delayEstimate+length(sweep)-1);
    
    % 计算相干性
    coherence_val = compute_mean_coherence(extracted, sweep, fs);
    
    if coherence_val < 0.6
        result.issues{end+1} = sprintf('信号相干性低(%.2f)<0.6', coherence_val);
    end
    
    % 计算SNR
    signal_energy = sum(extracted.^2);
    noise_est = median(abs(extracted - mean(extracted))) * 1.4826;
    snr_est = 20*log10(sqrt(signal_energy/length(extracted)) / (noise_est + eps));
    
    if snr_est < 10
        result.issues{end+1} = sprintf('信号SNR低(%.1fdB)<10dB', snr_est);
    end
    
    result.coherence = coherence_val;
    result.snr = snr_est;
else
    result.issues{end+1} = '延迟超出记录信号长度';
end

% ============== 综合判定 ==============
if isempty(result.issues)
    result.pass = true;
    result.confidence = 0.9;  % 高置信度
    validatedDelay = delayEstimate;
else
    % 存在问题，降低置信度
    numIssues = length(result.issues);
    result.confidence = max(0.1, 1 - numIssues * 0.2);
    
    % 如果失败但延迟在物理范围内，使用保守延迟
    if delayEstimate >= minPhysical && delayEstimate <= maxPhysical
        validatedDelay = delayEstimate;
    else
        % 使用物理中间值
        validatedDelay = round((minPhysical + maxPhysical) / 2);
        result.issues{end+1} = sprintf('使用保守延迟: %d', validatedDelay);
    end
end

result.delay_meters = distance_m;
result.validatedDelay = validatedDelay;
end