function metrics = assess_ir_quality(irData, snrData, peakPosData, coherenceEst, cfg)
% assess_ir_quality - 脉冲响应（IR）质量综合评估
%
% 输入:
%   irData        : [numSamples, numMics, numReps] 脉冲响应数据
%   snrData       : [numReps, numMics] 每次重复、每个麦克风的 SNR (dB)
%   peakPosData   : [numReps, numMics] 峰值位置（MATLAB 1-based 索引）
%   coherenceEst  : [numReps, numMics] 相干性估计值（0~1）
%   cfg           : 配置结构体（需包含 fs 等字段）
%
% 输出:
%   metrics       : 结构体，包含所有质量指标和可用于延迟估计的位置信息
%
% 注意:
%   - 峰值位置 peakPosData 必须是 1-based（MATLAB 默认）
%   - 延迟估计建议：delaySamples = (medianPeakPosPerMic - 1) - hardwareDelayCalibrated
%   - hardwareDelayCalibrated 应通过 loopback 校准获得（单位：samples）

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

% ===== 新增：每通道中值峰值位置（1-based index）=====
if numMics > 0 && ~isempty(peakPosMatrix)
    medianPeakPosPerMic = zeros(numMics, 1);
    for m = 1:numMics
        medianPeakPosPerMic(m) = median(peakPosMatrix(:, m));
    end
else
    medianPeakPosPerMic = zeros(numMics, 1);
end

% 峰值位置稳定性（全局）
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

% ==================== IR相似度分析（逐麦克风）====================
if numReps > 1
    minSimilarity = getCfgField(cfg, 'minSimilarity', 0.8);

    % 预分配
    nPairsPerMic = nchoosek(numReps, 2);
    totalPairs = numMics * nPairsPerMic;
    allMaxCorrs = zeros(totalPairs, 1);
    allZeroLagCorrs = zeros(totalPairs, 1);
    micAvgSim = zeros(numMics, 1);
    corrCache = zeros(numMics, numReps, numReps); % 上三角缓存
    
    idx = 1;
    for m = 1:numMics
        totalMaxCorr = 0;
        pairCount = 0;
        
        for i = 1:numReps
            for j = i+1:numReps
                ir1 = squeeze(irData(:, m, i));
                ir2 = squeeze(irData(:, m, j));
                
                [xc, ~] = xcorr(ir1, ir2, 'coeff');
                [maxCorr, ~] = max(xc);
                
                % 零延迟索引（MATLAB xcorr 输出长度 = 2*N-1，中心在 N）
                zeroLagIdx = floor(length(xc)/2) + 1;
                if zeroLagIdx <= length(xc)
                    zeroLagCorr = xc(zeroLagIdx);
                else
                    zeroLagCorr = 0;
                end
                
                allMaxCorrs(idx) = maxCorr;
                allZeroLagCorrs(idx) = zeroLagCorr;
                corrCache(m, i, j) = maxCorr;
                
                totalMaxCorr = totalMaxCorr + maxCorr;
                pairCount = pairCount + 1;
                idx = idx + 1;
            end
        end
        
        micAvgSim(m) = (pairCount > 0) ? totalMaxCorr / pairCount : 1;
    end
    
    avgSimilarity = mean(allMaxCorrs);
    avgZeroLagSimilarity = mean(allZeroLagCorrs);
    
    % 问题麦克风检测
    problematicMics = find(micAvgSim < (minSimilarity - 0.1));
    if ~isempty(problematicMics)
        fprintf('  警告: 麦克风 %s 相似度偏低 (<%.2f)，可能存在故障\n', ...
            mat2str(problematicMics'), minSimilarity - 0.1);
    end
    
    % 极性检查
    polarityIssues = sum(allMaxCorrs < -0.3);
    if polarityIssues > 0
        fprintf('  警告: 检测到%d个极性可能反转的测量对（相似度<%.1f）\n', ...
            polarityIssues, -0.3);
    end
    
    % 构建全局相似度矩阵
    similarityMatrix = zeros(numReps, numReps);
    for i = 1:numReps
        for j = i+1:numReps
            corrVals = corrCache(:, i, j);
            validCorr = corrVals(corrVals ~= 0);
            if ~isempty(validCorr)
                similarityMatrix(i, j) = mean(validCorr);
            else
                similarityMatrix(i, j) = 0;
            end
        end
    end
    similarityMatrix = similarityMatrix + similarityMatrix';
    similarityMatrix(1:numReps+1:end) = 1; % 对角线设为1

else
    minSimilarity = getCfgField(cfg, 'minSimilarity', 0.8);
    avgSimilarity = 1;
    avgZeroLagSimilarity = 1;
    similarityMatrix = 1;
    micAvgSim = ones(numMics, 1);
end
% ==================== 相似度分析结束 ====================

% 能量分布分析
irAvg = mean(irData, 3);
totalEnergy = sum(irAvg(:).^2) + eps;

% 前10%样本能量（预回声）
preSamples = max(1, floor(numSamples * 0.1));
preEnergy = sum(sum(irAvg(1:preSamples,:).^2));
preEnergyRatio = preEnergy / totalEnergy;

% 主能量窗口（±50 样本）
mainEnergy = 0;
for m = 1:numMics
    ir_mic = irAvg(:, m);
    [~, maxIdx] = max(abs(ir_mic));
    startIdx = max(1, maxIdx - 50);
    endIdx = min(numSamples, maxIdx + 50);
    mainEnergy = mainEnergy + sum(ir_mic(startIdx:endIdx).^2);
end
mainEnergyRatio = mainEnergy / totalEnergy;

% 相干性统计
if ~isempty(coherenceMatrix) && any(coherenceMatrix(:) > 0)
    meanCoherence = mean(coherenceMatrix(:));
    medianCoherence = median(coherenceMatrix(:));
    minCoherence = min(coherenceMatrix(:));
    
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

% 可用性判定阈值
snrThreshold = getCfgField(cfg, 'snrThresholdDB', 15);
coherenceThreshold = getCfgField(cfg, 'coherenceThreshold', 0.7);
maxPeakStd = getCfgField(cfg, 'maxPeakStd', 20);
maxPreEchoRatio = getCfgField(cfg, 'maxPreEchoRatio', 0.2);

% 处理单次测量
if numReps == 1
    stablePeaks = true;
else
    stablePeaks = (peakStd < maxPeakStd);
end

snrOK = (medianSNR >= snrThreshold);
lowPreEcho = (preEnergyRatio < maxPreEchoRatio);
consistent = (avgSimilarity > minSimilarity);
coherenceOK = (medianCoherence > coherenceThreshold);

% 综合评分
if numReps > 1
    weights = struct('snr', 0.3, 'coherence', 0.3, 'similarity', 0.2, ...
                     'stability', 0.1, 'preEcho', 0.1);
    score = weights.snr * double(snrOK) + ...
            weights.coherence * double(coherenceOK) + ...
            weights.similarity * double(consistent) + ...
            weights.stability * double(stablePeaks) + ...
            weights.preEcho * double(lowPreEcho);
    usable = (score >= 0.7) && snrOK && coherenceOK;
else
    usable = snrOK && coherenceOK && (preEnergyRatio < 0.1);
    score = double(usable);
end

% ==================== 输出结构 ====================
metrics = struct(...
    'medianSNR',             medianSNR, ...
    'peakStd',               peakStd, ...
    'peakRange',             peakRange, ...
    'stability',             stability, ...
    'similarity',            avgSimilarity, ...
    'zeroLagSimilarity',     avgZeroLagSimilarity, ...
    'preEnergyRatio',        preEnergyRatio, ...
    'mainEnergyRatio',       mainEnergyRatio, ...
    'meanCoherence',         meanCoherence, ...
    'medianCoherence',       medianCoherence, ...
    'minCoherence',          minCoherence, ...
    'snrOK',                 snrOK, ...
    'stablePeaks',           stablePeaks, ...
    'lowPreEcho',            lowPreEcho, ...
    'consistent',            consistent, ...
    'coherenceOK',           coherenceOK, ...
    'usable',                usable, ...
    'qualityScore',          score, ...
    'numSamples',            numSamples, ...
    'numMics',               numMics, ...
    'numReps',               numReps, ...
    'similarityMatrix',      similarityMatrix, ...
    'micAvgSimilarity',      micAvgSim, ...
    'medianPeakPosPerMic',   medianPeakPosPerMic);  % ✅ 关键新增：用于延迟估计

% ==================== 报告输出 ====================
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

fprintf('  3. 延迟相关信息:\n');
if numMics <= 5
    fprintf('     • 各通道中值峰值位置（1-based）: %s\n', mat2str(medianPeakPosPerMic'));
else
    fprintf('     • 各通道中值峰值位置（1-based）: [省略，共%d通道]\n', numMics);
end
fprintf('       → 物理延迟（样本）≈ (位置 - 1) - hardwareDelayCalibrated\n');

fprintf('  4. 综合评估:\n');
fprintf('     质量分数: %.2f/1.0\n', score);
fprintf('     判定结果: %s\n\n', ternary(usable, '✅ 测量质量合格', '❌ 测量质量不合格'));

if ~usable
    fprintf('  5. 改进建议:\n');
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

% 辅助函数：安全获取配置字段
function val = getCfgField(cfg, field, default)
    if isfield(cfg, field)
        val = cfg.(field);
    else
        val = default;
    end
end