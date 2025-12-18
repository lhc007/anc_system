function noise_rms_est = estimate_noise_floor(varargin)
% ESTIMATE_NOISE_FLOOR 实测麦克风本底噪声 RMS（用于 SNR 评估）
%
% 用法：
%   noise_rms_est = estimate_noise_floor(); 
%   noise_rms_est = estimate_noise_floor('Duration', 2.0, 'UseCache', true);
%
% 输入参数（可选 Name-Value）：
%   'Duration'    - 录制时长（秒），默认 2.0
%   'UseCache'    - 是否优先加载已保存的 noise_floor.mat，默认 false
%   'Verbose'     - 是否打印日志，默认 false
%
% 输出：
%   noise_rms_est - 误差麦克风通道噪声 RMS 的中位数（标量，归一化幅度）
%
% 依赖：
%   anc_config, hardware_init_measure

p = inputParser;
addParameter(p, 'Duration', 2.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'UseCache', false, @islogical);
addParameter(p, 'Verbose', false, @islogical);
parse(p, varargin{:});

durationSec = p.Results.Duration;
useCache    = p.Results.UseCache;
verbose     = p.Results.Verbose;

cfg = anc_config();
errMicIdx = cfg.micChannels.error;
fs = cfg.fs;
blockSize = cfg.timeFrameSamples;
totalSamples = round(durationSec * fs);
numBlocks = ceil(totalSamples / blockSize);

saveFile = 'calibration/noise_floor.mat';

% === 尝试从缓存加载 ===
if useCache && exist(saveFile, 'file')
    if verbose, fprintf('加载缓存噪声数据: %s\n', saveFile); end
    data = load(saveFile, 'noiseResult');
    if isfield(data, 'noiseResult') && ...
       isfield(data.noiseResult, 'noiseRmsEstRecommended')
        noise_rms_est = data.noiseResult.noiseRmsEstRecommended;
        return;
    end
end

% === 初始化硬件 ===
hw = [];
try
    if verbose, fprintf('初始化音频硬件（录制本底噪声）...\n'); end
    hw = hardware_init_measure(cfg);
    
    % 预热
    for pr = 1:cfg.preRollFrames
        hw.writer(zeros(blockSize, cfg.numSpeakers));
        hw.reader();
    end

    % === 录制静音 ===
    if verbose, fprintf('录制 %.1f 秒本底噪声...\n', durationSec); end
    recorded = zeros(numBlocks * blockSize, cfg.micNumChannels);
    wrPtr = 1;

    for b = 1:numBlocks
        hw.writer(zeros(blockSize, cfg.numSpeakers));
        micFrame = hw.reader();
        if isempty(micFrame)
            micFrame = zeros(blockSize, cfg.micNumChannels);
        elseif size(micFrame, 1) < blockSize
            micFrame(end+1:blockSize, :) = 0;
        end
        recorded(wrPtr:wrPtr+blockSize-1, :) = micFrame;
        wrPtr = wrPtr + blockSize;
    end
    recorded = recorded(1:totalSamples, :);

    % === 计算 RMS ===
    noiseRmsAll = zeros(cfg.micNumChannels, 1);
    for ch = 1:cfg.micNumChannels
        noiseRmsAll(ch) = rms(recorded(:, ch));
    end
    noiseRmsErr = noiseRmsAll(errMicIdx);
    noise_rms_est = median(noiseRmsErr);

    % === 保存结果 ===
    noiseResult = struct(...
        'fs', fs, ...
        'recordDurationSec', durationSec, ...
        'errorMicPhysicalChannels', errMicIdx, ...
        'noiseRmsErrorMics', noiseRmsErr, ...
        'noiseRmsEstRecommended', noise_rms_est, ...
        'timestamp', datestr(now) ...
    );
    mkdir(fileparts(saveFile));
    save(saveFile, 'noiseResult', '-v7.3');

    if verbose
        fprintf('✅ 噪声 RMS（误差麦中位数）: %.3e\n', noise_rms_est);
        fprintf('💾 已保存至: %s\n', saveFile);
    end

catch ME
    if verbose, fprintf('❌ 噪声测量失败: %s\n', ME.message); end
    rethrow(ME);
finally
    if ~isempty(hw)
        try, hw.release(); catch, end
    end
end
end