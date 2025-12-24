function test_speaker_channels()
% 测试双通道扬声器是否能独立发声
% 要求：听到左声道 → 右声道 → 左右同相 → 左右反相

clear; clc;

% === 配置（请根据你的实际 cfg 修改）===
cfg = anc_config();
fs = cfg.fs;
duration = 2.0;           % 每段 2 秒
f_tone = 800;             % 测试音频率
amp = 0.98;                % 幅度（避免削波）
deviceName = '扬声器2 (Realtek(R) Audio)';  % ← 替换为你的 cfg.spkDevice2Name

% === 创建 writer ===
w = audioDeviceWriter(...
    'Device', deviceName, ...
    'SampleRate', fs, ...
    'ChannelMappingSource', 'Property', ...
    'ChannelMapping', [1, 2]);  % 假设设备支持立体声

fprintf('▶ 准备测试扬声器通道...\n');
fprintf('  设备: %s\n', deviceName);
fprintf('  将依次播放：左声道 → 右声道 → 双声道同相 → 双声道反相\n');
pause(2);

t = (0:1/fs:(duration - 1/fs))';
tone = amp * sin(2*pi*f_tone*t);

% --- 1. 左声道 ---
fprintf('🔊 播放左声道 (通道1)...\n');
signal = [tone, zeros(size(tone))];  % [L, R]
play_full(w, signal, fs);

% --- 2. 右声道 ---
fprintf('🔊 播放右声道 (通道2)...\n');
signal = [zeros(size(tone)), tone];
play_full(w, signal, fs);

% --- 3. 双声道同相 ---
fprintf('🔊 播放双声道同相信号...\n');
signal = [tone, tone];
play_full(w, signal, fs);

% --- 4. 双声道反相（用于检测串扰）---
fprintf('🔊 播放双声道反相信号（理想情况下应减弱）...\n');
signal = [tone, -tone];
play_full(w, signal, fs);

release(w);
fprintf('✅ 测试完成。\n');

end

function play_full(writer, signal, fs)
% 将整个 signal 分块播放，避免内存溢出
blockSize = 1024;
numBlocks = ceil(size(signal,1) / blockSize);
for b = 1:numBlocks
    startIdx = (b-1)*blockSize + 1;
    endIdx = min(b*blockSize, size(signal,1));
    block = signal(startIdx:endIdx, :);
    % 补零到 blockSize（如果最后一块不足）
    if size(block,1) < blockSize
        block(end+1:blockSize, :) = 0;
    end
    step(writer, block);
end
pause(size(signal,1)/fs + 0.1);  % 等待播放完
end