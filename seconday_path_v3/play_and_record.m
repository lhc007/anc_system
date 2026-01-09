function [recorded, info] = play_and_record(hw, driveSig, spkIdx, cfg)
% PLAY_AND_RECORD - 同步播放单通道激励信号并录制多通道响应
%   driveSig: 单列向量 [N x 1]，要播放的激励信号（如 sweep）
%   spkIdx:   要播放到的扬声器索引（1-based）
%   recorded: 录制的麦克风信号 [M x numMics]

    totalSamples = length(driveSig);
    blockSize = cfg.timeFrameSamples;
    numBlocks = ceil(totalSamples / blockSize);

    % === 预热：清空音频缓冲区（防止前几帧污染）===
    for i = 1:cfg.preRollFrames
        hw.writer(zeros(blockSize, cfg.numSpeakers));
        hw.reader();
    end

    % === 分配足够大的录音缓冲区 ===
    tailSamples = round(2 * cfg.fs); % 多录2秒尾部残响
    maxRecordLen = numBlocks * blockSize + tailSamples + blockSize; % +1 block 安全余量
    recorded = zeros(maxRecordLen, cfg.micNumChannels);
    recPtr = 1;

    ptr = 1; % driveSig 指针

    % === 主播放+录音循环 ===
    for b = 1:numBlocks
        % 构造当前输出块（多通道）
        outBlock = zeros(blockSize, cfg.numSpeakers);
        if ptr <= totalSamples
            nCopy = min(blockSize, totalSamples - ptr + 1);
            outBlock(1:nCopy, spkIdx) = driveSig(ptr:ptr + nCopy - 1);
            ptr = ptr + nCopy;
        end

        % 播放并立即录音
        hw.writer(outBlock);
        micFrame = hw.reader();

        % --- 安全处理 micFrame 尺寸问题 ---
        frameRows = size(micFrame, 1);
        if frameRows == 0
            micFrame = zeros(blockSize, cfg.micNumChannels);
            frameRows = blockSize;
        end

        % 列数对齐
        if size(micFrame, 2) > cfg.micNumChannels
            micFrame = micFrame(:, 1:cfg.micNumChannels);
        elseif size(micFrame, 2) < cfg.micNumChannels
            micFrame = [micFrame, zeros(frameRows, cfg.micNumChannels - size(micFrame, 2))];
        end

        % 行数对齐（关键！某些声卡返回非 blockSize 行）
        if frameRows ~= blockSize
            if frameRows < blockSize
                % 补零到 blockSize
                micFrame = [micFrame; zeros(blockSize - frameRows, cfg.micNumChannels)];
            else
                % 截断（罕见，但安全起见）
                micFrame = micFrame(1:blockSize, :);
            end
        end

        % 存入 recorded
        recorded(recPtr:recPtr + blockSize - 1, :) = micFrame;
        recPtr = recPtr + blockSize;
    end

    % === 播放结束后继续录音（捕获残响）===
    tailBlocks = ceil(tailSamples / blockSize);
    for b = 1:tailBlocks
        hw.writer(zeros(blockSize, cfg.numSpeakers)); % 播放静音
        micFrame = hw.reader();

        frameRows = size(micFrame, 1);
        if frameRows == 0
            micFrame = zeros(blockSize, cfg.micNumChannels);
            frameRows = blockSize;
        end

        if size(micFrame, 2) > cfg.micNumChannels
            micFrame = micFrame(:, 1:cfg.micNumChannels);
        elseif size(micFrame, 2) < cfg.micNumChannels
            micFrame = [micFrame, zeros(frameRows, cfg.micNumChannels - size(micFrame, 2))];
        end

        if frameRows ~= blockSize
            if frameRows < blockSize
                micFrame = [micFrame; zeros(blockSize - frameRows, cfg.micNumChannels)];
            else
                micFrame = micFrame(1:blockSize, :);
            end
        end

        recorded(recPtr:recPtr + blockSize - 1, :) = micFrame;
        recPtr = recPtr + blockSize;
    end

    % === 裁剪未使用部分 ===
    recorded = recorded(1:recPtr - 1, :);

    % 返回信息
    info = struct(...
        'totalSamplesPlayed', totalSamples, ...
        'actualSamplesRecorded', size(recorded, 1), ...
        'blocksPlayed', numBlocks + tailBlocks);
end