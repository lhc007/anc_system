function metrics = assess_ir_quality(irData, snrData, peakPosData, coherenceEst, cfg)
% IR质量评估
% 计算多种质量指标
% 输入:
%   irData: [numSamples, numMics, numReps] 脉冲响应数据
%   snrData: [numReps, numMics] SNR估计值
%   peakPosData: [numReps, numMics] 峰值位置
%   coherenceEst: [numReps, numMics] 相干性估计值
%   cfg: 配置结构体
% 输出:
%   metrics: 质量指标结构体

fprintf('[assess_ir_quality] 开始IR质量评估\n');
fprintf('  数据维度: %d样本 × %d麦克风 × %d重复\n', ...
    size(irData,1), size(irData,2), size(irData,3));
fprintf('  相干性矩阵大小: %dx%d\n', size(coherenceEst,1), size(coherenceEst,2));

% 获取数据维度
[numSamples, numMics, numReps] = size(irData);

% 基本统计
snrMatrix = snrData;
peakPosMatrix = peakPosData;
coherenceMatrix = coherenceEst;

% 中值SNR
medianSNR = median(snrMatrix(:));

% 峰值位置稳定性
if numReps > 1 && numMics > 1
    peakStd = std(peakPosMatrix(:));
    peakRange = max(peakPosMatrix(:)) - min(peakPosMatrix(:));
else
    peakStd = 0;
    peakRange = 0;
end

% 标准化稳定性指标
if numReps > 1
    stability = 1 - min(1, peakStd / 100); % 基于峰值标准差
else
    stability = 0.95; % 单次测量，假设稳定性较高但不完美
end

% ==================== 改进：IR相似度（重复间一致性）====================
% 使用归一化互相关（NCC）的最大值作为相似度指标
if numReps > 1
    similarityMatrix = zeros(numReps, numReps);
    similarityAtZeroLag = zeros(numReps, numReps); 
    validPairs = 0;
    
    for i = 1:numReps
        for j = i+1:numReps
            % 对所有麦克风取平均，得到单通道代表IR
            ir1 = mean(squeeze(irData(:, :, i)), 2);  % [N x 1]
            ir2 = mean(squeeze(irData(:, :, j)), 2);  % [N x 1]
            
            % 确保长度一致（安全处理）
            minLen = min(length(ir1), length(ir2));
            ir1 = ir1(1:minLen);
            ir2 = ir2(1:minLen);
            
            % 计算归一化互相关 (Normalized Cross-Correlation)
            [xc, lags] = xcorr(ir1, ir2, 'coeff');
            
            % 找到最大互相关值（即最佳对齐下的相似度）
            [maxCorr, ~] = max(xc);

            % 获取零延迟处的互相关（检查时间对齐）- 更高效的方法
            zeroLagIdx = round(length(xc)/2) + 1;  % xcorr输出的中间位置是零延迟
            if zeroLagIdx <= length(xc)
                zeroLagCorr = xc(zeroLagIdx);
            else
                zeroLagCorr = 0;
            end
            
            % 存储相似度
            similarityMatrix(i, j) = maxCorr;
            similarityAtZeroLag(i, j) = zeroLagCorr;
            validPairs = validPairs + 1;
        end
    end
    
    if validPairs > 0
        avgSimilarity = sum(similarityMatrix(:)) / validPairs;
        avgZeroLagSimilarity = sum(similarityAtZeroLag(:)) / validPairs;
        
        % 警告：如果最佳对齐相似度远高于零延迟相似度，说明存在时间偏移
        if (avgSimilarity - avgZeroLagSimilarity) > 0.2
            fprintf('  警告: 测量间存在明显时间偏移(%.3f vs %.3f)\n', ...
                avgSimilarity, avgZeroLagSimilarity);
        end
    else
        avgSimilarity = 0;
        avgZeroLagSimilarity = 0;
    end

    % ==================== 检查极性一致性 ====================
    % 排除对角线元素，只检查下三角
    mask = tril(true(size(similarityMatrix)), -1);
    similarityValues = similarityMatrix(mask);
    polarityIssues = sum(similarityValues < -0.3); % -0.3阈值更敏感
    
    if polarityIssues > 0
        fprintf('  警告: 检测到%d个极性可能反转的测量对（相似度<%.1f）\n', ...
            polarityIssues, -0.3);
    end
    
    % 对称化矩阵以便存储
    similarityMatrix = similarityMatrix + similarityMatrix';
    for k = 1:numReps
        similarityMatrix(k, k) = 1;
    end
else
    avgSimilarity = 1;  % 单次测量，相似度为1
    avgZeroLagSimilarity = 1;
    similarityMatrix = 1;  % 1x1矩阵
end
% ==================== 相似度计算结束 ====================

% 能量分布
irAvg = mean(irData, 3);
totalEnergy = sum(irAvg(:).^2);

% 前10%样本的能量占比（检查预回声）
preSamples = max(1, floor(numSamples * 0.1));
preEnergy = sum(sum(irAvg(1:preSamples,:).^2));
preEnergyRatio = preEnergy / (totalEnergy + 1e-12);

% 主能量窗口（找到峰值窗口）
if numMics > 0
    % 找到所有IR的最大绝对值位置
    maxPos = zeros(numMics, 1);
    
    for m = 1:numMics
        ir_mic = mean(irData(:, m, :), 3);
        [~, maxPos(m)] = max(abs(ir_mic));
    end
    
    % 主能量窗口（峰值周围±50个样本）
    mainEnergy = 0;
    for m = 1:numMics
        startIdx = max(1, maxPos(m) - 50);
        endIdx = min(numSamples, maxPos(m) + 50);
        mainEnergy = mainEnergy + sum(irAvg(startIdx:endIdx, m).^2);
    end
    
    mainEnergyRatio = mainEnergy / (totalEnergy + 1e-12);
else
    mainEnergyRatio = 0;
end

% 相干性统计
if ~isempty(coherenceMatrix) && any(coherenceMatrix(:) > 0)
    meanCoherence = mean(coherenceMatrix(:));
    medianCoherence = median(coherenceMatrix(:));
    minCoherence = min(coherenceMatrix(:));
    
    % 检查相干性是否异常高（接近1.0）
    if medianCoherence > 0.99
        fprintf('  警告: 相干性异常高(%.3f)，可能存在计算错误\n', medianCoherence);
    elseif medianCoherence < 0.3
        fprintf('  警告: 相干性过低(%.3f)，测量质量可能不佳\n', medianCoherence);
    end
else
    meanCoherence = 0;
    medianCoherence = 0;
    minCoherence = 0;
end

% 可用性判定
% 从cfg获取阈值，如果不存在则使用默认值
if isfield(cfg, 'snrThresholdDB')
    snrThreshold = cfg.snrThresholdDB;
else
    snrThreshold = 15; % 默认15dB
end

if isfield(cfg, 'coherenceThreshold')
    coherenceThreshold = cfg.coherenceThreshold;
else
    coherenceThreshold = 0.7; % 默认0.7
end

if isfield(cfg, 'maxPeakStd')
    maxPeakStd = cfg.maxPeakStd;
else
    maxPeakStd = 20; % 默认20样本
end

if isfield(cfg, 'minSimilarity')
    minSimilarity = cfg.minSimilarity;
else
    minSimilarity = 0.8; % 默认0.8
end

if isfield(cfg, 'maxPreEchoRatio')
    maxPreEchoRatio = cfg.maxPreEchoRatio;
else
    maxPreEchoRatio = 0.2; % 预回声能量占比小于20%
end

% 处理单次测量时的特殊情况
if numReps == 1
    % 单次测量时，相似度阈值应该更宽松
    minSimilarity = 0.5;
    % 峰值稳定性不适用于单次测量
    stablePeaks = true;  % 设为true，因为无法评估
else
    stablePeaks = (peakStd < maxPeakStd);
end

snrOK = (medianSNR >= snrThreshold);
lowPreEcho = (preEnergyRatio < maxPreEchoRatio);
consistent = (avgSimilarity > minSimilarity);
coherenceOK = (medianCoherence > coherenceThreshold);

% ==================== 改进：更精细的可用性判定 ====================
% 根据重复次数调整权重
if numReps > 1
    % 多次测量：为不同指标分配权重
    weights = struct(...
        'snr', 0.3, ...
        'coherence', 0.3, ...
        'similarity', 0.2, ...
        'stability', 0.1, ...
        'preEcho', 0.1);
    
    score = 0;
    score = score + weights.snr * double(snrOK);
    score = score + weights.coherence * double(coherenceOK);
    score = score + weights.similarity * double(consistent);
    score = score + weights.stability * double(stablePeaks);
    score = score + weights.preEcho * double(lowPreEcho);
    
    % 阈值可以根据实际应用调整
    if score >= 0.7 && snrOK && coherenceOK
        usable = true;
    else
        usable = false;
    end
else
    % 单次测量：更严格的要求
    usable = snrOK && coherenceOK && (preEnergyRatio < 0.1); % 单次时要求更低的预回声
    score = double(usable);  % 单次测量时质量分数就是0或1
end

% ==================== 构建输出结构（扩展）====================
metrics = struct(...
    'medianSNR', medianSNR, ...
    'peakStd', peakStd, ...
    'peakRange', peakRange, ...
    'stability', stability, ...
    'similarity', avgSimilarity, ...
    'zeroLagSimilarity', avgZeroLagSimilarity, ... % 新增
    'preEnergyRatio', preEnergyRatio, ...
    'mainEnergyRatio', mainEnergyRatio, ...
    'meanCoherence', meanCoherence, ...
    'medianCoherence', medianCoherence, ...
    'minCoherence', minCoherence, ...
    'snrOK', snrOK, ...
    'stablePeaks', stablePeaks, ...
    'lowPreEcho', lowPreEcho, ...
    'consistent', consistent, ...
    'coherenceOK', coherenceOK, ...
    'usable', usable, ...
    'qualityScore', score, ... % 新增：量化质量分数
    'numSamples', numSamples, ...
    'numMics', numMics, ...
    'numReps', numReps, ...
    'similarityMatrix', similarityMatrix); % 新增：完整的相似度矩阵

% ==================== 改进输出报告 ====================
fprintf('\n  ===== IR质量评估报告 =====\n');
fprintf('  1. 基本参数:\n');
fprintf('     样本数: %d, 麦克风: %d, 重复: %d\n', numSamples, numMics, numReps);

fprintf('  2. 质量指标:\n');
fprintf('     • SNR: %.1f dB %s\n', medianSNR, ternary(snrOK, '(合格)', '(不合格)'));
fprintf('     • 相干性: %.3f %s\n', medianCoherence, ternary(coherenceOK, '(合格)', '(不合格)'));
fprintf('     • 相似度: %.3f (最佳对齐) / %.3f (零延迟) %s\n', ...
    avgSimilarity, avgZeroLagSimilarity, ternary(consistent, '(合格)', '(不合格)'));

if numReps > 1
    fprintf('     • 峰值稳定性: σ=%.1f样本, 范围=%d样本 %s\n', ...
        peakStd, peakRange, ternary(stablePeaks, '(合格)', '(不合格)'));
else
    fprintf('     • 峰值稳定性: 单次测量无法评估 %s\n', ...
        ternary(stablePeaks, '(合格)', '(不合格)'));
end

fprintf('     • 能量分布: 主能量=%.1f%%, 预回声=%.1f%% %s\n', ...
    mainEnergyRatio*100, preEnergyRatio*100, ternary(lowPreEcho, '(合格)', '(不合格)'));

fprintf('  3. 综合评估:\n');
fprintf('     质量分数: %.2f/1.0\n', score);
fprintf('     判定结果: %s\n\n', ternary(usable, '✅ 测量质量合格', '❌ 测量质量不合格'));

if ~usable
    fprintf('  4. 改进建议:\n');
    if ~snrOK
        fprintf('     • 提高信噪比: 增大激励信号幅度或降低环境噪声\n');
    end
    if ~coherenceOK
        fprintf('     • 改善相干性: 检查非线性失真、系统同步问题\n');
    end
    if ~consistent
        fprintf('     • 提高一致性: 确保测量条件稳定，检查时间对齐\n');
    end
    if numReps > 1 && ~stablePeaks
        fprintf('     • 提高稳定性: 减少环境扰动，固定测量装置\n');
    end
    if ~lowPreEcho
        fprintf('     • 减少预回声: 检查系统延迟，优化反卷积算法\n');
    end
end
end
