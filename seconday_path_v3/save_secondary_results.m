function save_secondary_results(impulseResponses, cfg, exciteInfo, allSNRs, allCohs)
    fs = cfg.fs;
    errMicIdx = cfg.micChannels.error;
    Lh = cfg.irMaxLen;
    numErrMics = length(errMicIdx);
    
    secondary = struct();
    secondary.impulseResponses = impulseResponses;
    secondary.fs = fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;
    secondary.errorMicPhysicalChannels = errMicIdx;
    secondary.irLength = Lh;
    secondary.version = 'v5.1';
    secondary.timestamp = datetime('now', 'TimeZone', 'local');
    secondary.measurementDate = datetime('now', 'TimeZone', 'local');
    
    % 添加元数��
    secondary.config = struct();
    fieldList = {'fs', 'numSpeakers', 'irMaxLen', 'repetitions', ...
        'snrThresholdDB', 'coherenceThreshold', 'sweepFreqStartHz', ...
        'sweepFreqEndHz', 'deconvPreDelayKeep'};
    
    for f = fieldList
        if isfield(cfg, f{1})
            secondary.config.(f{1}) = cfg.(f{1});
        end
    end
    
    secondary.excitationInfo = exciteInfo;
    secondary.measurementSNRs = allSNRs;
    secondary.measurementCoherences = allCohs;

    % 处理目标路径与文件名，保证是字符向量
    [saveDir, name, ext] = fileparts(cfg.secondaryPathFile);
    if isempty(ext)
        ext = '.mat';
    end
    if isempty(saveDir)
        saveDir = pwd; % 若未指定目录，使用当前工作目录
    end
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    saveFilePath = fullfile(saveDir, [name ext]);  % char

    % 若文件存在则添加时间戳（使用安全的 datetime -> char）
    if exist(saveFilePath, 'file')
        ts = char(datetime('now', 'TimeZone', 'local', 'Format', 'yyyyMMdd_HHmmss'));
        newFilename = fullfile(saveDir, [name '_' ts ext]);
        fprintf('[warning] 文件已存在，保存为: %s\n', newFilename);
        saveFilePath = newFilename;
    end

    % 最终保存（确保传入的是 char）
    saveFilePath = char(saveFilePath);
    save(saveFilePath, 'secondary', '-v7.3');
    fprintf('\n[success] 次级路径已保存至: %s\n', saveFilePath);
    
    % 显示测量总结
    show_measurement_summary_final(impulseResponses, allSNRs, allCohs, cfg);
end

function show_measurement_summary_final(impulseResponses, allSNRs, allCohs, cfg)
    fprintf('\n=== 测量质量总结 ===\n');
    
    numSpk = size(impulseResponses, 3);
    numMics = size(impulseResponses, 2);
    
    total_measurements = numSpk * numMics * cfg.repetitions;
    % 计算时忽略 NaN
    good_snr = sum(~isnan(allSNRs(:)) & allSNRs(:) >= cfg.snrThresholdDB);
    good_coh = sum(~isnan(allCohs(:)) & allCohs(:) >= cfg.coherenceThreshold);
    
    fprintf('总测量次数: %d\n', total_measurements);
    fprintf('SNR合格率: %d/%d (%.1f%%)\n', good_snr, total_measurements, ...
        good_snr/total_measurements*100);
    fprintf('相干性合格率: %d/%d (%.1f%%)\n', good_coh, total_measurements, ...
        good_coh/total_measurements*100);
    
    for spk = 1:numSpk
        fprintf('\n扬声器 %d:\n', spk);
        for mic = 1:numMics
            snr_vals = squeeze(allSNRs(:, mic, spk));
            coh_vals = squeeze(allCohs(:, mic, spk));
            % 使用 nanmedian/nanstd 忽略无效重复
            fprintf('  麦克风 %d: SNR=%.1f±%.1f dB, Coh=%.3f±%.3f (基于有效重复)\n', ...
                mic, median(snr_vals), std(snr_vals), ...
                median(coh_vals), std(coh_vals));
        end
    end
end