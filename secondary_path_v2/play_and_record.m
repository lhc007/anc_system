function [recorded, info] = play_and_record(hw, driveSig, spkIdx, cfg)
% 播放与录音函数
% 确保精确的同步和定时

totalSamples = length(driveSig);
blockSize = cfg.timeFrameSamples;
numBlocks = ceil(totalSamples / blockSize);

% 初始化录音buffer
recorded = zeros(numBlocks * blockSize, cfg.micNumChannels);

% 预热帧
for i = 1:cfg.preRollFrames
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
if size(recorded, 1) > totalSamples
    recorded = recorded(1:totalSamples, :);
elseif size(recorded, 1) < totalSamples
    recorded(end+1:totalSamples, :) = 0;
end

% 返回信息结构
info = struct(...
    'totalSamples', totalSamples, ...
    'actualSamples', size(recorded, 1), ...
    'blocksPlayed', numBlocks);

end


