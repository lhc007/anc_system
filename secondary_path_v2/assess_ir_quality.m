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
    avgSimilarity = 1;  % 单次测量，相似度为1
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
    snrThreshold = 15; % 默认15dB（提高要求）
end

if isfield(cfg, 'coherenceThreshold')
    coherenceThreshold = cfg.coherenceThreshold;
else
    coherenceThreshold = 0.7; % 默认0.7
end

if isfield(cfg, 'maxPeakStd')
    maxPeakStd = cfg.maxPeakStd;
else
    maxPeakStd = 20; % 默认20样本（更严格）
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

snrOK = (medianSNR >= snrThreshold);
stablePeaks = (peakStd < maxPeakStd);
lowPreEcho = (preEnergyRatio < maxPreEchoRatio);
consistent = (avgSimilarity > minSimilarity);
coherenceOK = (medianCoherence > coherenceThreshold);

% 综合可用性判定
if numReps > 1
    % 多次测量：允许部分指标不通过
    % 使用加权判定
    criteria_count = 5;  % 总共5个标准
    passed_count = sum([snrOK, stablePeaks, lowPreEcho, consistent, coherenceOK]);
    pass_ratio = passed_count / criteria_count;
    
    % 如果通过率>60%且关键指标（SNR和相干性）通过，则可用
    usable = (pass_ratio > 0.6) && snrOK && coherenceOK;
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


%% 相干性低通常意味着：
 
% 环境噪声大（风扇、空调、电脑底噪）
% 非线性失真（扬声器过载、削波）
% 系统不稳定（播放/录制不同步、缓冲区溢出）
% 扫频激励能量不足（幅度过小，SNR 虽够但相干性仍低）
% 反卷积或对齐有误差

%% 输出详细质量报告
fprintf('  质量评估结果:\n');
fprintf('    - SNR: %.1f dB (阈值: %.1f dB) %s\n', ...
    medianSNR, snrThreshold, ternary(snrOK, '✓', '✗'));
fprintf('    - 相干性: %.2f (阈值: %.2f) %s\n', ...
    medianCoherence, coherenceThreshold, ternary(coherenceOK, '✓', '✗'));
fprintf('    - 峰值稳定性: %.1f 样本 (阈值: <%.1f) %s\n', ...
    peakStd, maxPeakStd, ternary(stablePeaks, '✓', '✗'));
fprintf('    - 预回声比例: %.3f (阈值: <%.2f) %s\n', ...
    preEnergyRatio, maxPreEchoRatio, ternary(lowPreEcho, '✓', '✗'));
fprintf('    - 相似度: %.2f (阈值: >%.2f) %s\n', ...
    avgSimilarity, minSimilarity, ternary(consistent, '✓', '✗'));
fprintf('    - 综合判定: %s\n', ternary(usable, '可用', '不可用'));
end


function str = ternary(condition, trueStr, falseStr)
% 三元操作符模拟
if condition
    str = trueStr;
else
    str = falseStr;
end
end


