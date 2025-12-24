% record_test_mic_advanced.m
% 录制六通道麦克风，并仅保存通道 5 和 6
clear; clc;

%% 配置
fs           = 48000;
durationSec  = 8;
deviceName   = '六通道麦克风阵列 (YDM6MIC Audio)';
targetChannels = [5, 6];        % ←← 指定要保留的通道（1-based）
nSaveChannels  = length(targetChannels);  % = 2

filenameOut = 'recorded_ch5_ch6.wav';

% 必须读取全部 6 通道（因为设备从 CH1 开始输出）
reader = audioDeviceReader( ...
    'Device',         deviceName, ...
    'SampleRate',     fs, ...
    'NumChannels',    6, ...      % ← 必须是 6
    'SamplesPerFrame',512 ...
);

totalSamples = fs * durationSec;
% 预分配：只存目标通道（2列）
recorded = zeros(totalSamples, nSaveChannels);

% 预热
for k = 1:5; reader(); end

fprintf('🎙️ 准备录制 %d 秒（保存通道 %s）...\n', durationSec, mat2str(targetChannels));
pause(1);
disp('🔴 开始录音！请立即播放 sweep 音频。');

samplesSoFar = 0;
while samplesSoFar < totalSamples
    frameAll = reader();  % [N x 6]
    available = min(size(frameAll,1), totalSamples - samplesSoFar);
    % 只提取 CH5 和 CH6
    recorded(samplesSoFar+1 : samplesSoFar+available, :) = frameAll(1:available, targetChannels);
    samplesSoFar = samplesSoFar + available;
end

release(reader);

% 保存双通道 WAV
audiowrite(filenameOut, recorded, fs);
fprintf('✅ 录音完成: %s（%d 通道）\n', filenameOut, nSaveChannels);