function hw = hardware_init_measure(cfg)
% hardware_init_measure
% 初始化用于次级路径测量的音频硬件（麦克风 + 扬声器）
%
% 输出 hw 结构体包含：
%   - .reader()    : 无参函数，返回 [SamplesPerFrame x numInputs] 音频帧
%   - .writer(block): 接收 [N x numOutputs] 音频块进行播放
%   - .release()   : 释放设备
%   - .numInputs   : 麦克风通道数
%   - .fs          : 采样率
%   - .stats       : 错误统计

fprintf('[hardware-measure] 初始化硬件...\n');

% 确保必要字段存在
if ~isfield(cfg, 'timeFrameSamples'), cfg.timeFrameSamples = 160; end
if ~isfield(cfg, 'micNumChannels'),   cfg.micNumChannels = 1; end
if ~isfield(cfg, 'fs'),               cfg.fs = 48000; end

hw = struct();
hw.stats = struct('readErrors', 0, 'writeErrors', 0);
hw.numInputs = cfg.micNumChannels;  % ←←← 关键补充！供 play_and_record 使用
hw.fs = cfg.fs;                     % ←←← 便于其他函数获取采样率

% === 麦克风读取器 ===
try
    rdr = audioDeviceReader( ...
        'Device', cfg.micDeviceName, ...
        'NumChannels', cfg.micNumChannels, ...
        'SampleRate', cfg.fs, ...
        'SamplesPerFrame', cfg.timeFrameSamples);  
catch ME
    error('[hardware-measure] audioDeviceReader 创建失败: %s', ME.message);
end

% === 扬声器写入器 ===
% 注意：ChannelMapping 根据实际扬声器数量配置
try
    w = audioDeviceWriter(...
        'Device', cfg.spkDevice2Name, ...
        'SampleRate', cfg.fs, ...
        'ChannelMappingSource', 'Property', ...
        'ChannelMapping', 1:cfg.numSpeakers); % 映射到物理输出通道 1～2
catch ME
    try release(rdr); catch; end
    error('[hardware-measure] audioDeviceWriter 创建失败: %s\n', ME.message);
end

% 封装接口（使用函数句柄隐藏底层对象）
hw.reader  = @() safeRead(rdr, cfg, hw);
hw.writer  = @(block) safeWrite(w, block, hw);
hw.release = @() releaseAll(rdr, w);

fprintf('[hardware-measure] 设备初始化完成:\n');
fprintf('  麦克风设备: %s (通道数: %d)\n', cfg.micDeviceName, cfg.micNumChannels);
fprintf('  扬声器设备: %s\n', cfg.spkDevice1Name);
fprintf('  采样率: %d Hz, 帧长: %d 样本\n', cfg.fs, cfg.timeFrameSamples);
end

% --- 辅助函数 ---
function frame = safeRead(rdr, cfg, hw)
    try
        % 兼容新旧 MATLAB 版本：R2016b+ 支持直接调用对象
        if verLessThan('matlab', '9.1') % R2016b 之前用 step()
            frame = step(rdr);
        else
            frame = rdr();
        end
    catch
        hw.stats.readErrors = hw.stats.readErrors + 1;
        frame = zeros(cfg.timeFrameSamples, cfg.micNumChannels, 'like', 0);
    end
    % 确保输出维度正确（防止驱动返回不足通道）
    if size(frame, 2) < cfg.micNumChannels
        frame = [frame, zeros(size(frame,1), cfg.micNumChannels - size(frame,2), 'like', frame)];
    elseif size(frame, 2) > cfg.micNumChannels
        frame = frame(:, 1:cfg.micNumChannels);
    end
end

function safeWrite(wObj, block, hw)
    try
        if verLessThan('matlab', '9.1')
            step(wObj, block);
        else
            wObj(block);
        end
    catch
        hw.stats.writeErrors = hw.stats.writeErrors + 1;
    end
end

function releaseAll(rdr, w)
    try 
        if verLessThan('matlab', '9.1')
            release(rdr); release(w);
        else
            rdr.release(); w.release();
        end
    catch
        % 忽略释放错误
    end
    fprintf('[hardware-measure] 音频设备已释放\n');
end