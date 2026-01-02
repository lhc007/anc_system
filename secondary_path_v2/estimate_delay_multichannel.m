function [allDelays, allConfidences, allCorrelations] = estimate_delay_multichannel(irAligned, sweepCore_scaled, cfg, numErrMics)
% 多通道延迟估计
fs = cfg.fs();

% 使用平均对齐信号
allRecorded = zeros(size(irAligned{1}));
for rep = 1:length(irAligned)
    allRecorded = allRecorded + irAligned{rep};
end
avgRecorded = allRecorded / length(irAligned);

% 对每个麦克风通道估计延迟
allDelays = zeros(numErrMics, 1);
allConfidences = zeros(numErrMics, 1);
allCorrelations = zeros(numErrMics, 1);

for m = 1:numErrMics
    fprintf('  通道 %d 延迟估计... ', m);
    
    [delayEst, valResult, ~] = robust_delay_estimation_pipeline(...
        avgRecorded(:, m), sweepCore_scaled, cfg);
    
    allDelays(m) = delayEst;
    allConfidences(m) = valResult.confidence;
    allCorrelations(m) = valResult.correlation;
    
    fprintf('延迟: %d 样本 (置信度: %.2f)\n', delayEst, valResult.confidence);
end
end


function [delayEstimate, validationResult, diagnostics] = robust_delay_estimation_pipeline(...
    avgRecorded_aligned, sweepCore, cfg)
% 鲁棒延迟估计管道（适配click对齐后）
% avgRecorded_aligned: 已对齐的录制信号（长度=sweepCore）
% sweepCore: 核心扫频信号（参考）
% sweepCore_ref: 同上（冗余，保持接口兼容）
% cfg: 配置
% spk_info: 扬声器信息（可选）

fs = cfg.fs;
diagnostics = struct();

%% 由于已对齐，延迟估计简化为峰值检测
% 计算互相关（使用对齐后的信号）
[corr_vals, lags] = xcorr(avgRecorded_aligned, sweepCore, 'normalized');
positive_lags = lags >= 0;
if any(positive_lags)
    corr_vals_pos = corr_vals(positive_lags);
    lags_pos = lags(positive_lags);
    
    [max_corr, max_idx] = max(corr_vals_pos);
    peak_lag = lags_pos(max_idx);
    
    % 峰值质量评估
    noise_region = corr_vals_pos;
    peak_region = max(max_idx-10, 1):min(max_idx+10, length(corr_vals_pos));
    noise_region(peak_region) = [];
    
    if ~isempty(noise_region)
        noise_level = median(abs(noise_region)) * 1.4826; % MAD
    else
        noise_level = 0.01;
    end
    
    peak_snr = max_corr / (noise_level + eps);
    peak_quality = max_corr * sqrt(peak_snr); % 综合指标
    
    fprintf('对齐后延迟估计: lag=%d, 峰值=%.3f, SNR=%.1f\n', peak_lag, max_corr, 20*log10(peak_snr));
else
    peak_lag = 0;
    max_corr = 0;
    peak_quality = 0;
    fprintf('警告: 未找到有效峰值\n');
end

%% 延迟验证
% 使用峰值质量评估
if peak_quality > 0.5 && max_corr > 0.3
    validation_pass = true;
    confidence = min(0.9, peak_quality);
elseif peak_quality > 0.2 && max_corr > 0.1
    validation_pass = true;
    confidence = peak_quality;
else
    validation_pass = false;
    confidence = 0.2;
end

%% 构建结果
delayEstimate = peak_lag; % 对齐后，延迟通常为0或很小
validationResult = struct(...
    'pass', validation_pass, ...
    'confidence', confidence, ...
    'delay_global', delayEstimate, ...
    'correlation', max_corr, ...
    'snr_db', 20*log10(peak_snr + eps));
diagnostics.peak_lag = peak_lag;
diagnostics.peak_correlation = max_corr;
diagnostics.peak_snr = peak_snr;
diagnostics.peak_quality = peak_quality;

fprintf('最终延迟估计: %d 样本 (置信度: %.1f%%)\n', delayEstimate, confidence*100);
end


