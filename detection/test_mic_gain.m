%% 实时显示麦克风输入幅度
% test_mic_gain.m
fs = 48000;
reader = audioDeviceReader('Device', '六通道麦克风阵列 (YDM6MIC Audio)', ...
    'SampleRate', fs, 'NumChannels', 6, 'SamplesPerFrame', 1024);

fprintf('正在监听麦克风... 请对着麦克风说话或播放声音\n');
for i = 1:100
    x = reader();
    rms_val = rms(x,1);
    fprintf('通道 RMS: ');
    fprintf('%.2e ', rms_val);
    fprintf('\n');
    pause(0.1);
end
close(reader);