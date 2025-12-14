function hw = hardware_init_measure(cfg)
% hardware_init_measure
% 初始测量硬件

if ~isfield(cfg,'timeFrameSamples'), cfg.timeFrameSamples = 160; end

hw = struct();
hw.stats = struct('readErrors',0,'writeErrors',0);

% === 麦克风读取器（通常用 WASAPI，保持默认）===
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
try
    w = audioDeviceWriter('Device', cfg.spkDevice2Name,'SampleRate', cfg.fs);
catch ME
    try release(rdr); catch; end
    error('[hardware-measure] audioDeviceWriter (ASIO) 创建失败: %s\n可能原因：\n1. ASIO4ALL 未安装\n2. 设备未在 ASIO4ALL 中启用\n3. 其他程序占用了音频设备', ME.message);
end

% 封装接口
hw.reader  = @() safeRead(rdr, cfg, hw);
hw.writer  = @(block) safeWrite(w, block, hw);   % ← 单一 writer
hw.release = @() releaseAll(rdr, w);

fprintf('[hardware-measure] 设备初始化完成:\n');
fprintf('  麦克风: %s\n', cfg.micDeviceName);
fprintf('  扬声器: %s (ASIO, 4ch)\n', cfg.spkDevice2Name);
end

% --- 辅助函数 ---
function frame = safeRead(rdr, cfg, hw)
    try
        frame = step(rdr);
    catch
        hw.stats.readErrors = hw.stats.readErrors + 1;
        frame = zeros(cfg.timeFrameSamples, cfg.micNumChannels);
    end
    if size(frame,2) < cfg.micNumChannels
        frame(:, end+1:cfg.micNumChannels) = 0;
    end
end

function safeWrite(wObj, block, hw)
    try
        step(wObj, block);
    catch
        hw.stats.writeErrors = hw.stats.writeErrors + 1;
    end
end

function releaseAll(rdr, w)
    try release(rdr); catch; end
    try release(w); catch; end
    fprintf('[hardware-measure] 设备已释放\n');
end