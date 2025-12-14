%% test_speaker_for_anc.m
% 测试扬声器播放强度
clear; clc;

fs = 48000;
duration = 1.0;        % 每个频率播放 1 秒
testFreqs = [100, 80, 60, 50, 40];  % 从高到低测试
amplitude = 0.85;      % 足够大但不削波

% 使用系统默认播放设备（最可靠）
player = audioDeviceWriter('SampleRate', fs);

fprintf('🔊 测试扬声器输出能力（ANC 次级路径）...\n');
fprintf('请靠近扬声器仔细听每个频率的声音。\n\n');

for k = 1:length(testFreqs)
    f = testFreqs(k);
    t = (0:1/fs:duration)';
    signal = amplitude * sin(2*pi*f*t);  % 列向量！
    
    fprintf('▶ 播放 %d Hz ... ', f);
    player(signal);
    pause(duration + 0.2);
    fprintf('完成\n');
end

release(player);
fprintf('\n✅ 测试结束。如果 60 Hz 以下几乎无声，则不适用于低频 ANC。\n');