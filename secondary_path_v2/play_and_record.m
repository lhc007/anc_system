function [recorded, info] = play_and_record(hw, driveSig, spkIdx, cfg)
% 播放与录音函数
% 确保精确的同步和定时

totalSamples = length(driveSig);
blockSize = cfg.timeFrameSamples;
numBlocks = ceil(totalSamples / blockSize);

totalSamples = length(driveSig);
redundancy = round(1.5 * cfg.fs);
recorded = zeros(numBlocks * blockSize + redundancy, cfg.micNumChannels);

% 预热帧
for i = 1:cfg.preRollFrames
    hw.writer(zeros(blockSize, cfg.numSpeakers));
    hw.reader();
end

% === 新增：额外静音， ===
for i = 1:2
    hw.writer(zeros(blockSize, cfg.numSpeakers));
    hw.reader();
end

% 主循环
ptr = 1;
recPtr = 1;

for b = 1:numBlocks
    % 准备输出块
    outBlock = zeros(blockSize, cfg.numSpeakers);
    
    % 填充当前块的激励信号
    if ptr <= totalSamples
        nToCopy = min(blockSize, totalSamples - ptr + 1);
        outBlock(1:nToCopy, spkIdx) = driveSig(ptr:ptr+nToCopy-1);
        ptr = ptr + nToCopy;
    end
    
    % 播放并录音
    hw.writer(outBlock);
    micFrame = hw.reader();
    
    % 处理可能的尺寸不匹配
        if isempty(micFrame)
            micFrame = zeros(blockSize, cfg.micNumChannels);
        elseif size(micFrame, 2) ~= cfg.micNumChannels
            % 如果列数不匹配，截断或填充
            if size(micFrame, 2) > cfg.micNumChannels
                micFrame = micFrame(:, 1:cfg.micNumChannels);
            else
                % 填充到正确列数
                micFrame = [micFrame, zeros(blockSize, cfg.micNumChannels - size(micFrame, 2))];
            end
        end
    
    % 保存录音
    recorded(recPtr:recPtr+blockSize-1, :) = micFrame;
    recPtr = recPtr + blockSize;
end

% 截断到正确长度
actualRecordedLength = recPtr - 1; % 已写入的样本数
if actualRecordedLength < size(recorded, 1)
    recorded = recorded(1:actualRecordedLength, :); % 只裁掉未使用的尾部零
end

% 返回实际录制长度
info = struct(...
    'totalSamplesPlayed', totalSamples, ...
    'actualSamplesRecorded', size(recorded, 1), ...
    'blocksPlayed', numBlocks);
end


