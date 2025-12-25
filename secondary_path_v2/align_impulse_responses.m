function irAligned = align_impulse_responses(irData, peakPosData, cfg)
% IR对齐函数
% 基于峰值位置对齐所有重复

[numSamples, numMics, numReps] = size(irData);

% 计算参考峰值位置
allPeaks = peakPosData(:);
if ~isempty(allPeaks)
    refPeak = round(median(allPeaks));
else
    refPeak = cfg.minPhysDelaySamples;
end

% 对齐每个重复
irAligned = zeros(size(irData));

for rep = 1:numReps
    % 计算该重复的平均峰值
    repPeak = round(mean(peakPosData(rep, :)));
    
    % 计算偏移量
    shift = refPeak - repPeak;
    
    for m = 1:numMics
        irOrig = squeeze(irData(:, m, rep));
        
        if shift > 0
            % 向右移动（在左侧补零）
            irShifted = [zeros(shift, 1); irOrig(1:end-shift)];
        elseif shift < 0
            % 向左移动（在右侧补零）
            irShifted = [irOrig(-shift+1:end); zeros(-shift, 1)];
        else
            irShifted = irOrig;
        end
        
        % 确保长度一致
        if length(irShifted) > numSamples
            irShifted = irShifted(1:numSamples);
        elseif length(irShifted) < numSamples
            irShifted = [irShifted; zeros(numSamples - length(irShifted), 1)];
        end
        
        irAligned(:, m, rep) = irShifted;
    end
end
end