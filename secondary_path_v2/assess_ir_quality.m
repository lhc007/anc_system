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
    peakRange = range(peakPosMatrix(:));
else
    peakStd = 0;
    peakRange = 0;
end

% 标准化稳定性指标
if numReps > 1
    stability = 1 - min(1, peakStd / 100); % 基于峰值标准差
else
    stability = 1; % 单次测量，稳定性为1
end

% IR相似度（重复间一致性）
if numReps > 1
    similarityMatrix = zeros(numReps, numReps);
    validPairs = 0;
    
    for i = 1:numReps
        for j = i+1:numReps
            % 计算所有麦克风的平均IR
            ir1 = mean(squeeze(irData(:,:,i)), 2);
            ir2 = mean(squeeze(irData(:,:,j)), 2);
            
            % 归一化
            ir1_norm = ir1 - mean(ir1);
            ir2_norm = ir2 - mean(ir2);
            
            % 计算相关系数
            numerator = sum(ir1_norm .* ir2_norm);
            denominator = sqrt(sum(ir1_norm.^2) * sum(ir2_norm.^2));
            
            if denominator > 1e-12
                corrVal = numerator / denominator;
                similarityMatrix(i,j) = corrVal;
                validPairs = validPairs + 1;
            else
                similarityMatrix(i,j) = 0;
            end
        end
    end
    
    if validPairs > 0
        avgSimilarity = sum(similarityMatrix(:)) / validPairs;
    else
        avgSimilarity = 0;
    end
else
    avgSimilarity = 1; % 单次测量，相似度为1
end

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
    maxVals = zeros(numMics, 1);
    maxPos = zeros(numMics, 1);
    
    for m = 1:numMics
        ir_mic = mean(irData(:, m, :), 3);
        [maxVals(m), maxPos(m)] = max(abs(ir_mic));
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
if ~isempty(coherenceMatrix)
    meanCoherence = mean(coherenceMatrix(:));
    medianCoherence = median(coherenceMatrix(:));
    minCoherence = min(coherenceMatrix(:));
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
    snrThreshold = 10; % 默认10dB
end

if isfield(cfg, 'coherenceThreshold')
    coherenceThreshold = cfg.coherenceThreshold;
else
    coherenceThreshold = 0.7; % 默认0.7
end

if isfield(cfg, 'maxPeakStd')
    maxPeakStd = cfg.maxPeakStd;
else
    maxPeakStd = 50; % 默认50样本
end

if isfield(cfg, 'minSimilarity')
    minSimilarity = cfg.minSimilarity;
else
    minSimilarity = 0.8; % 默认0.8
end

snrOK = (medianSNR >= snrThreshold);
stablePeaks = (peakStd < maxPeakStd);
lowPreEcho = (preEnergyRatio < 0.2); % 预回声能量占比小于20%
consistent = (avgSimilarity > minSimilarity);
coherenceOK = (medianCoherence > coherenceThreshold);

% 综合可用性判定
if numReps > 1
    % 多次测量：要求所有条件都满足
    usable = snrOK && stablePeaks && lowPreEcho && consistent && coherenceOK;
else
    % 单次测量：只要求SNR和相干性
    usable = snrOK && coherenceOK;
end

% 构建输出结构
metrics = struct(...
    'medianSNR', medianSNR, ...
    'peakStd', peakStd, ...
    'peakRange', peakRange, ...
    'stability', stability, ...
    'similarity', avgSimilarity, ...
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
    'numSamples', numSamples, ...
    'numMics', numMics, ...
    'numReps', numReps);
end