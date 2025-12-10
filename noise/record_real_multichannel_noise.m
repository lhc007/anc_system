%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% record_real_multichannel_noise.m — 纯录制模式
% 录制初级噪声场（即没有 ANC 控制信号时的原始噪声传播）
% 使用外部设备（如手机）播放 road_noise.wav，
% 本机仅用6通道麦克风录制真实的初级噪声场。
% 输出：anc_sim_data_road_noise.wav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function record_real_multichannel_noise()
    cfg = anc_config();
    
    % === 检查录制时长 ===
    if ~isfield(cfg, 'duration_sec') || cfg.duration_sec < 1
        error('❌ cfg.duration_sec 必须 ≥ 1 秒（建议 10~30 秒）');
    end
    
    total_samples = round(cfg.duration_sec * cfg.fs);
    num_frames = ceil(total_samples / cfg.frameSize);
    total_samples = num_frames * cfg.frameSize; % 对齐帧边界

    % === 尝试打开设备 ===
    try
        recorder = audioDeviceReader(...
            'Device', cfg.micDeviceName, ...
            'NumChannels', 6, ...
            'SampleRate', cfg.fs, ...
            'SamplesPerFrame', cfg.frameSize);
    catch ME
        error('❌ 无法打开录音设备 "%s"：\n%s', cfg.micDeviceName, ME.message);
    end

    % 初始化缓冲区
    recorded = zeros(total_samples, 6);
    fprintf('🎙️ 开始纯录制（请确保外部设备正在播放 noise）...\n');
    fprintf(' 录制设备: %s\n', cfg.micDeviceName);
    fprintf(' 采样率: %d Hz, 通道数: %d\n', recorder.SampleRate, recorder.NumChannels);
    fprintf(' 录制时长: %.1f 秒 (%d 样本)\n', cfg.duration_sec, total_samples);
    fprintf(' ⏳ 请在 %d 秒内保持外部噪声播放！\n', cfg.duration_sec);

    % 清空初始不稳定数据
    for i = 1:10, recorder(); end

    % 主录制循环
    try
        for f = 0:num_frames-1
            audioIn = recorder(); % [frameSize x 6]
            startIdx = f * cfg.frameSize + 1;
            recorded(startIdx:startIdx+cfg.frameSize-1, :) = audioIn;
            if mod(f, max(1, floor(num_frames/10))) == 0
                fprintf(' 进度: %.0f%%\n', 100*(f+1)/num_frames);
            end
        end
    catch
        fprintf('\n🛑 录制被中断。\n');
    end

    clear recorder; % 释放设备

    % === 检查信号能量 ===
    if max(abs(recorded(:))) < 1e-6
        warning('⚠️ 警告：录制信号能量极低！请检查麦克风硬件。');
    else
        fprintf('🔊 录制信号峰值: %.2e\n', max(abs(recorded(:))));
    end

    % 保存文件
    audiowrite(cfg.outputFile, recorded, cfg.fs);
    fprintf('✅ 录制完成！文件已保存为: %s\n', cfg.outputFile);
end