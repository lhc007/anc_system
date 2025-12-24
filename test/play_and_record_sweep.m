% play_and_record_sweep.m
% ✅ 修复版：确保左→右顺序播放 sweep，并完整录制 CH5 & CH6
% 兼容 MATLAB R2025a
% 播放设备：扬声器2 (Realtek(R) Audio)
% 录制设备：六通道麦克风阵列 (YDM6MIC Audio)

clear; clc;
cfg = anc_config();
%% 配置参数
fs            = cfg.fs;                     % 采样率
sweepFile     = 'test_sweep_48k_5s.wav';   % 输入 sweep 文件
filenameOut   = 'recorded_ch5_ch6.wav';    % 输出录音文件

playbackDevice = cfg.spkDevice2Name;
recordDevice   = cfg.micDeviceName;

targetChannels = cfg.micChannels.error;                   % 只录 CH5 和 CH6
nSaveChannels  = length(targetChannels);

sweepCfg = struct('fs',cfg.fs,'T',cfg.sweepDuration,'f1',cfg.sweepF1,'f2',cfg.sweepF2,...
    'padLeading',cfg.padLeading,'padTrailing',cfg.padTrailing,'amplitude',cfg.amplitude);
% %% 加载 sweep 信号
% [sweep, fs_sweep] = audioread(sweepFile);
% if fs_sweep ~= fs
%     error('Sweep 采样率 (%d) 与配置采样率 (%d) 不匹配！', fs_sweep, fs);
% end
% 
% sweep = sweep(:);              % 确保为列向量


%% === 1. 生成 sweep ===
fprintf('🔄 正在生成 ESS 扫频信号...\n');
[sweepSig, ~, sweepInfo] = generate_sweep(sweepCfg);
if sweepInfo.fs ~= fs
    error('生成的 sweep 采样率 (%d) 与系统配置采样率 (%d) 不匹配！', sweepInfo.fs, fs);
end
sweep = sweepSig(:); % 确保为列向量
fprintf('✅ 生成完成: %.1f–%.0f Hz, %.2f 秒\n', sweepInfo.f1, sweepInfo.f2, sweepInfo.fullLength/fs);

sweepLen = length(sweep);
gapSamples = round(0.5 * fs);  % 0.5 秒静音（样本数）

%% === 2. 构造 playSignal (左 → 静音 → 右) ===
% 左声道播放 sweep，右=0
leftBlock  = [sweep, zeros(sweepLen, 1)];
% 中间静音
silence    = zeros(gapSamples, 2);
% 右声道播放 sweep，左=0
rightBlock = [zeros(sweepLen, 1), sweep];

% 拼接完整播放序列
playSignal = [leftBlock; silence; rightBlock];  % 总长度: sweepLen + gap + sweepLen
totalPlaySamples = size(playSignal, 1);
totalPlayTimeSec = totalPlaySamples / fs;

% 录制总时长：播放时间 + 3秒余量（防截断）
recordDurationSec = totalPlayTimeSec + 3;
totalRecordSamples = round(recordDurationSec * fs);

fprintf('🔊 播放总时长: %.2f 秒\n', totalPlayTimeSec);
fprintf('🎙️ 录制总时长: %.2f 秒\n', recordDurationSec);

%% 初始化音频设备
fprintf('⚙️ 正在初始化音频硬件...\n');
hw = hardware_init_measure(cfg);

player = hw.writer;   % 函数句柄：@(block) safeWrite(...)
reader = hw.reader;   % 函数句柄：@() safeRead(...)


%% 预热设备（发送几帧静音）
warmupFrames = 5;
for k = 1:warmupFrames
    silentFrame = zeros(512, 2);
    player(silentFrame);
    reader();
end

fprintf('\n✅ 设备预热完成。\n');
pause(1);
disp('🔴 开始播放（左 → 静音 → 右）并同步录制...');

%% 同步播放 + 录制
recorded = zeros(totalRecordSamples, nSaveChannels);
samplesRecorded = 0;

% 按帧播放完整 playSignal
frameSize = 512;
totalPlayFrames = ceil(totalPlaySamples / frameSize);

for f = 1:totalPlayFrames
    startIdx = (f - 1) * frameSize + 1;
    endIdx   = min(f * frameSize, totalPlaySamples);
    playFrame = playSignal(startIdx:endIdx, :);
    
    % 补零至 frameSize（最后一帧可能不足）
    if size(playFrame, 1) < frameSize
        playFrame(end+1:frameSize, :) = 0;
    end
    
    player(playFrame);
    
    % 立即读取录制数据
    recFrameAll = reader();  % 512 x 6
    nRec = size(recFrameAll, 1);
    
    % 存入录制缓冲区
    if samplesRecorded + nRec <= totalRecordSamples
        recorded(samplesRecorded+1 : samplesRecorded+nRec, :) = ...
            recFrameAll(:, targetChannels);
        samplesRecorded = samplesRecorded + nRec;
    else
        % 最后一截
        remaining = totalRecordSamples - samplesRecorded;
        if remaining > 0
            recorded(samplesRecorded+1 : end, :) = ...
                recFrameAll(1:remaining, targetChannels);
            samplesRecorded = totalRecordSamples;
        end
    end
end

% 继续录制剩余时间（确保录满）
while samplesRecorded < totalRecordSamples
    recFrameAll = reader();
    nRec = size(recFrameAll, 1);
    remaining = totalRecordSamples - samplesRecorded;
    if remaining <= 0
        break;
    end
    useRec = min(nRec, remaining);
    recorded(samplesRecorded+1 : samplesRecorded+useRec, :) = ...
        recFrameAll(1:useRec, targetChannels);
    samplesRecorded = samplesRecorded + useRec;
end

%% 清理资源
try hw.release(); catch; end

%% 保存录制结果
audiowrite(filenameOut, recorded, fs);
fprintf('\n✅ 播放与录制完成！\n');
fprintf('📁 录音已保存为: %s\n', filenameOut);
fprintf('📊 录制样本数: %d (%.2f 秒)\n', size(recorded,1), size(recorded,1)/fs);