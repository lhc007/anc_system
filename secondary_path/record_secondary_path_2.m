% record_secondary_path.m 多扬声器→多误差麦克风次级路径测量（ESS）
% 优化点：多麦对齐、IR峰值对齐、内存控制、动态SNR、错误安全释放
% 版本：v2.6 优化增强版（2025-12-22）

clear; clc;

%% ============== 配置区 ==============
cfg = anc_config();

% 仅处理误差麦克风[5 6]
errMicIdx = cfg.micChannels.error;
numErrMics = length(errMicIdx);
fprintf('[measure] 误差麦克风通道: %s\n', mat2str(errMicIdx));

% Sweep 配置
sweepCfg = struct('fs',cfg.fs,'T',cfg.sweepDuration,'f1',cfg.sweepF1,'f2',cfg.sweepF2,...
    'padLeading',cfg.padLeading,'padTrailing',cfg.padTrailing,'amplitude',cfg.amplitude);

% deconvolve 参数
deconvParams = struct( ...
    'extraTail',     cfg.deconvExtraTail, ...
    'preDelayKeep',  cfg.deconvPreDelayKeep, ...
    'peakThreshDB',  cfg.deconvPeakThreshDB, ...
    'tailTotal',     cfg.irMaxLen, ...
    'maxSearch',     cfg.deconvMaxSearch , ...
    'regEps',        cfg.deconvRegEps, ...
    'noiseWin',      cfg.deconvNoiseWin, ...
    'envSmoothWin',  cfg.deconvEnvSmoothWin , ...
    'cumEnergyFrac', cfg.deconvCumEnergyFrac, ...
    'minPeakFrac',   cfg.deconvMinPeakFrac, ...
    'snrBodyRadius', cfg.deconvSnrBodyRadius, ...
    'fftCorrEnable', cfg.deconvFftCorrEnable, ...
    'debugMode',     cfg.deconvDebugMode, ...
    'peakRefineEnable', true, ...
    'peakRefineRadius', cfg.peakRefineRadius, ...
    'minPhysDelay',  cfg.minPhysDelaySamples,...  
    'maxPhysDelay',  cfg.maxPhysDelaySamples,  ... 
    'delaySearchRadius', cfg.delaySearchRadius ... 
    );

%% ============== 初始化硬件（带安全释放）=============
hw = []; % 初始化为空，便于 finally 安全释放
try
    %% ============== 生成 Sweep ==============
    fprintf('[measure] 生成 ESS sweep...\n');
    [sweepSigCore, invFilter, sweepInfo] = generate_sweep(sweepCfg);  % 无 padding 的核心 sweep

    % 显式保存 core（用于后续 xcorr 和 deconv）
    sweepCore = sweepSigCore(:);  % 列向量
    sweepCoreLen = length(sweepCore);

    % 构造带 silence 的完整播放信号
    preSilence  = zeros(round(cfg.preSilenceSec * cfg.fs), 1);
    postSilence = zeros(round(cfg.postSilenceSec * cfg.fs), 1);

    if cfg.enableLowFreqBoost
        [b_lp,a_lp] = butter(4, cfg.lowFreqCutHz/(cfg.fs/2), 'low');
        sweepLow    = filter(b_lp,a_lp,sweepCore);
        sweepMixed  = (1 - cfg.lowFreqMixRatio)*sweepCore + cfg.lowFreqMixRatio*sweepLow;
        sweepMixed  = sweepMixed / (max(abs(sweepMixed)) + 1e-12) * sweepCfg.amplitude;
        sweepDrive  = sweepMixed;
    else
        sweepDrive  = sweepCore / (max(abs(sweepCore)) + 1e-12) * sweepCfg.amplitude;
    end

    sweepSig = [preSilence; sweepDrive; postSilence];
    totalSamples = length(sweepSig);

    % 记录 core 在完整信号中的位置（可选，用于调试）
    sweepCoreStart = length(preSilence) + 1;
    sweepCoreEnd   = sweepCoreStart + sweepCoreLen - 1;
    fprintf('[measure] Sweep 总长=%d 样本, Sweep Core 长度=%d\n', totalSamples, sweepCoreLen);

    %% ============== 容器初始化 ==============
    Lh = cfg.irMaxLen;
    impulseResponses = zeros(Lh, numErrMics, cfg.numSpeakers);
    impulseResponsesAligned = zeros(Lh, numErrMics, cfg.numSpeakers);

    blockSize = cfg.timeFrameSamples;
    numBlocks = ceil(totalSamples / blockSize);

    % 默认不存 raw，除非明确开启（节省内存）
    saveRaw = cfg.saveAllRaw || cfg.saveFirstRaw;
    if saveRaw
        rawRecordings = cell(cfg.numSpeakers, cfg.repetitions);
    else
        rawRecordings = [];
    end
    
    % 保存每个repeat的IR用于诊断
    if cfg.saveEachRepeatIR
        irEachRepeat = cell(cfg.numSpeakers, 1);
    end
    
    meta = struct();
    meta.global = struct('sampleRate',cfg.fs,'sweepLen',totalSamples,'repetitions',cfg.repetitions,...
        'minPhysDelaySamples',cfg.minPhysDelaySamples,'maxPhysDelaySamples',cfg.maxPhysDelaySamples,...
        'enableRepeatAlignment',cfg.enableRepeatAlignment,'lowFreqBoost',cfg.enableLowFreqBoost,...
        'exportAlignedIR',cfg.exportAlignedIR);
    meta.perSpeaker = cell(cfg.numSpeakers,1);
    allIrRepList = cell(cfg.numSpeakers,1);

    %% ============== 测量循环 ==============
    for spk = 1:cfg.numSpeakers
        fprintf('\n[measure] ===== 扬声器 %d =====\n', spk);
        repChirpStarts = zeros(cfg.repetitions, 1);

        repGlobalShifts = zeros(cfg.repetitions,1);
        repRecordedRaw = cell(cfg.repetitions,1);
        repDiscardRepeat = false(cfg.repetitions,1);
        repPeakLists_raw = cell(cfg.repetitions,1);
        repPeakLists_forIR = cell(cfg.repetitions,1);
        repPeakReliable = cell(cfg.repetitions,1);
        repSNRlist = zeros(cfg.repetitions, numErrMics);
        repEnergyConcentration = zeros(cfg.repetitions, numErrMics);
        repPeakSharpness = zeros(cfg.repetitions, numErrMics);
        irRepList = zeros(Lh, numErrMics, cfg.repetitions);

        % === 阶段1：采集所有 repeats ===
        for rep = 1:cfg.repetitions
            fprintf('[measure]  播放重复 %d/%d...\n', rep, cfg.repetitions);
            if rep == 1
                hw = hardware_init_measure(cfg);
            else
                % 释放上一次的硬件
                if ~isempty(hw) && ismethod(hw, 'release')
                    hw.release();
                end
                pause(0.3);  % 给 Windows 音频子系统时间清理资源
                hw = hardware_init_measure(cfg);  % 重新打开
            end

            % === 连续指针模式播放 ===
            spkDriveSig = cfg.spkAmplitude(spk) * sweepSig(:);  % 列向量，全长
            
            % 防削波处理
            maxVal = max(abs(spkDriveSig));
            if maxVal > 0.98
                scaleFactor = 0.95 / maxVal;
                spkDriveSig = spkDriveSig * scaleFactor;
                fprintf('[measure] 播放缩放 %.2f\n', scaleFactor);
            else
                scaleFactor = 1.0;  % 显式定义，即使未缩放
            end
            
            coreStartIdx = length(preSilence) + 1;
            actualPlayedCore = spkDriveSig(coreStartIdx : coreStartIdx + sweepCoreLen - 1);
            refSweepForXcorr = actualPlayedCore;
            
            % 初始化指针和录音 buffer
            ptr = 1;
            recordedFull = zeros(numBlocks * blockSize, cfg.micNumChannels);

            % 预热
            for pr = 1:cfg.preRollFrames
                hw.writer(zeros(blockSize, cfg.numSpeakers));
                hw.reader();
            end

            % 主播放循环（连续指针）
            rdPtr = 1;
            for b = 1:numBlocks
                outBlock = zeros(blockSize, cfg.numSpeakers);
                if ptr <= totalSamples
                    nToCopy = min(blockSize, totalSamples - ptr + 1);
                    outBlock(1:nToCopy, spk) = spkDriveSig(ptr : ptr + nToCopy - 1);
                    ptr = ptr + nToCopy;
                end
                % 播放并录制
                hw.writer(outBlock);
                micFrame = hw.reader();
                if isempty(micFrame)
                    micFrame = zeros(blockSize, cfg.micNumChannels);
                elseif size(micFrame, 1) < blockSize
                    micFrame(end+1:blockSize, :) = 0;
                end
                recordedFull(rdPtr : rdPtr + blockSize - 1, :) = micFrame;
                rdPtr = rdPtr + blockSize;
            end

            % 截取与原始 sweep 等长（确保后续互相关对齐）
            recordedFull = recordedFull(1:totalSamples, :);
            recorded = recordedFull(:, errMicIdx);

             % 实时监控信号质量
            if cfg.enableRealTimeMonitor
                monitor_measurement_quality(recorded, spkDriveSig, cfg.fs, spk, rep);
            end

            if saveRaw && ((spk==1 && rep==1 && cfg.saveFirstRaw) || cfg.saveAllRaw)
                rawRecordings{spk,rep} = recorded;
            end
            repRecordedRaw{rep} = recorded;

             % === 互相关估计延迟 ===
             [chirpStartEst, correlationQuality] = estimate_delay_robust_multichannel(...
                recorded, refSweepForXcorr, cfg.fs, ...
                cfg.minPhysDelaySamples, cfg.maxPhysDelaySamples);

             % 在互相关估计后添加调试信息
            fprintf('[DEBUG] recorded size: %d x %d\n', size(recorded));
            fprintf('[DEBUG] refSweepForXcorr size: %d\n', length(refSweepForXcorr));
            fprintf('[DEBUG] recorded RMS: %.4f\n', rms(recorded(:)));
            fprintf('[DEBUG] refSweep RMS: %.4f\n', rms(refSweepForXcorr));

            if ~isfinite(chirpStartEst) || isnan(chirpStartEst) || correlationQuality < 0.1
                chirpStartEst = round((cfg.minPhysDelaySamples + cfg.maxPhysDelaySamples) / 2);
                fprintf('[WARN] Spk%d Rep%d: 延迟估计失败 → 使用默认偏移 = %d\n', ...
                    spk, rep, chirpStartEst);
                correlationQuality = 0;
            end

            repChirpStarts(rep) = chirpStartEst;
            fprintf('[measure]  Spk%d Rep%d chirpStart=%d (corr=%.3f)\n', ...
                spk, rep, chirpStartEst, correlationQuality);
            
            if correlationQuality < 0.3
                fprintf('[WARN]  Spk%d Rep%d 互相关质量低 (%.3f < 0.3)\n', ...
                    spk, rep, correlationQuality);
            end
        end

        % === 阶段2：确定参考 chirp 起始点（用于诊断，非窗口中心）===
        validShifts = repChirpStarts(~isnan(repChirpStarts) & isfinite(repChirpStarts));
        if ~isempty(validShifts)
            if length(validShifts) > 2
                validShiftsSorted = sort(validShifts);
                trimmed = validShiftsSorted(2:end-1);
                refChirpStart = round(mean(trimmed));
            else
                refChirpStart = round(median(validShifts));
            end
        else
            refChirpStart = cfg.minPhysDelaySamples + sweepCoreLen/2;
        end
        fprintf('[diag] Spk%d 参考 chirp 起始位置: %d 样本\n', spk, refChirpStart);

       % === 阶段3：窗口截取 + 反卷积 ===
        for rep = 1:cfg.repetitions
            if repDiscardRepeat(rep)
                continue;
            end

            rec = repRecordedRaw{rep};
            chirpStartEst = repChirpStarts(rep);

            if isnan(chirpStartEst) || ~isfinite(chirpStartEst)
                chirpStartEst = refChirpStart;
            end

            % === 关键：安全截取窗口（chirp 起始对齐）===
            [recExtract, offset, windowInfo] = extract_window_safely(rec, chirpStartEst, ...
                sweepCoreLen, cfg.irMaxLen);

            % === 关键修复：绝不补零！若 chirp 不完整，跳过 ===
            if size(recExtract, 1) < sweepCoreLen
                fprintf('[WARN] Spk%d Rep%d: chirp not fully captured (len=%d < %d). Skipping.\n', ...
                    spk, rep, size(recExtract,1), sweepCoreLen);
                repDiscardRepeat(rep) = true;
                % 初始化无效 IR
                irCurrent = zeros(Lh, numErrMics);
                peakIdxEach_raw = NaN(numErrMics, 1);
                peakRelEach = false(numErrMics, 1);
                snrMic = -Inf(numErrMics, 1);
                energyConcentration = NaN(numErrMics, 1);
                peakSharpness = NaN(numErrMics, 1);
            else
                if windowInfo.warningFlag
                    fprintf('[WARN] Spk%d Rep%d: %s\n', spk, rep, windowInfo.warningMsg);
                end

                sweepForDeconv = sweepCore;
                irCurrent = zeros(Lh, numErrMics);
                peakIdxEach_raw = zeros(numErrMics, 1);
                peakRelEach = false(numErrMics, 1);
                snrMic = zeros(numErrMics, 1);
                energyConcentration = zeros(numErrMics, 1);
                peakSharpness = zeros(numErrMics, 1);
                irQuality = cell(numErrMics, 1);

                for m = 1:numErrMics
                    outStruct = deconvolve_sweep_v2(recExtract(:,m), sweepForDeconv, cfg.fs, deconvParams);
                    h_full = outStruct.h;
                    pk_local_global = outStruct.pkLocalGlobal;
                    pk_global = offset + pk_local_global;

                    if isfield(outStruct, 'irQuality')
                        irQuality{m} = outStruct.irQuality;
                    end

                    % 能量集中度 & 峰值锐度
                    if ~isempty(h_full) && length(h_full) > 10
                        energySeq = abs(h_full).^2;
                        totalE = sum(energySeq);
                        if totalE > 1e-18
                            cumE = cumsum(energySeq) / totalE;
                            energy85 = find(cumE >= 0.85, 1, 'first');
                            energyConcentration(m) = energy85;
                            
                            ir_abs = abs(h_full);
                            peakVal = max(ir_abs);
                            validStart = max(1, min(length(h_full), cfg.minPhysDelaySamples+1));
                            validEnd = min(length(h_full), validStart + Lh - 1);
                            if validEnd >= validStart
                                validIR = h_full(validStart:validEnd);
                                rmsValid = sqrt(mean(abs(validIR).^2));
                                if rmsValid > 0
                                    peakSharpness(m) = peakVal / rmsValid;
                                end
                            end
                        end
                    end

                    % 若主峰太早，尝试用能量累积重定位
                    if pk_local_global < cfg.minPhysDelaySamples && ~isempty(h_full)
                        energySeq = abs(h_full).^2;
                        totalE = sum(energySeq);
                        if totalE > 1e-18
                            cumE = cumsum(energySeq) / totalE;
                            newPk_local = find(cumE >= cfg.deconvCumEnergyFrac, 1, 'first');
                            if ~isempty(newPk_local) && newPk_local > pk_local_global
                                pk_local_global = newPk_local;
                                pk_global = offset + pk_local_global;
                            end
                        end
                    end

                    inDelayRange = (pk_local_global >= cfg.minPhysDelaySamples) && ...
                                   (pk_local_global <= cfg.maxPhysDelaySamples);
                    hasGoodSNR = outStruct.snrEst >= cfg.snrThresholdDB;
                    hasReliablePeak = isfield(outStruct, 'peakReliability') && ...
                                     outStruct.peakReliability;
                    
                    reliable = inDelayRange && hasGoodSNR && hasReliablePeak;
                    
                    if cfg.deconvDebugMode
                        delay_status = 'FAIL'; if inDelayRange, delay_status = 'OK'; end
                        snr_status   = 'FAIL'; if hasGoodSNR,   snr_status = 'OK'; end
                        fprintf('  [DEBUG] Mic%d: pk_global=%d (范围[%d,%d] → %s), ', ...
                            m, pk_global, cfg.minPhysDelaySamples, cfg.maxPhysDelaySamples, delay_status);
                        fprintf('SNR=%.2f dB (thresh=%.1f → %s), ', outStruct.snrEst, cfg.snrThresholdDB, snr_status);
                        fprintf('peakRel=%d → reliable=%d\n', hasReliablePeak, reliable);
                    end

                    peakIdxEach_raw(m) = pk_global;
                    peakRelEach(m) = reliable;
                    snrMic(m) = outStruct.snrEst;
                    
                    Luse = min(Lh, length(h_full));
                    irCurrent(1:Luse, m) = h_full(1:Luse);
                end
            end

            repPeakLists_raw{rep} = peakIdxEach_raw;
            repPeakReliable{rep} = peakRelEach;
            repSNRlist(rep, :) = snrMic;
            repEnergyConcentration(rep, :) = energyConcentration;
            repPeakSharpness(rep, :) = peakSharpness;
            irRepList(:, :, rep) = irCurrent;
        end
        
        if cfg.saveEachRepeatIR
            irEachRepeat{spk} = irRepList;
        end

         % === 4. 基于 IR 峰值做重复间对齐（若启用）===
        allReliablePeaks = [];
        for rep = 1:cfg.repetitions
            if repDiscardRepeat(rep), continue; end
            peaks = repPeakLists_raw{rep}(repPeakReliable{rep});
            allReliablePeaks = [allReliablePeaks; peaks];
        end

       if isempty(allReliablePeaks)
            refPeakForAlign = cfg.minPhysDelaySamples + 50;
            fprintf('[WARN] Spk%d 没有可靠的峰值，使用默认对齐参考: %d\n', spk, refPeakForAlign);
        else
            refPeakForAlign = round(median(allReliablePeaks));
            fprintf('[diag] Spk%d 对齐参考峰值: %d (基于 %d 个可靠峰值)\n', ...
                spk, refPeakForAlign, length(allReliablePeaks));
        end

        % 统一 IR 集合用于平均
        if cfg.doAlignRepeats || cfg.exportAlignedIR
            irToAverage = zeros(Lh, numErrMics, cfg.repetitions);
            for rep = 1:cfg.repetitions
                if repDiscardRepeat(rep)
                    irToAverage(:, :, rep) = zeros(Lh, numErrMics);
                    continue;
                end
                ir_orig = irRepList(:, :, rep);
                reliableMask = repPeakReliable{rep};
                if any(reliableMask)
                    repPeak = mean(repPeakLists_raw{rep}(reliableMask));
                else
                    repPeak = refPeakForAlign;
                end
                shiftIR = round(refPeakForAlign - repPeak);
                maxShift = min(100, Lh/4);
                if abs(shiftIR) > maxShift
                    fprintf('[WARN] Spk%d Rep%d 对齐偏移过大: %d，限制为 %d\n', ...
                        spk, rep, shiftIR, sign(shiftIR)*maxShift);
                    shiftIR = sign(shiftIR) * maxShift;
                end

                for m = 1:numErrMics
                    h = ir_orig(:, m);
                    if shiftIR > 0
                        h_shifted = [zeros(shiftIR,1); h(1:end-shiftIR)];
                    elseif shiftIR < 0
                        h_shifted = [h(-shiftIR+1:end); zeros(-shiftIR,1)];
                    else
                        h_shifted = h;
                    end
                    irToAverage(:, m, rep) = h_shifted(1:Lh);
                end
            end
        else
            irToAverage = irRepList;
        end

        % === 5. 平均 IR（跳过丢弃的重复）===
        validReps = ~repDiscardRepeat;
        if any(validReps)
            irAvg = mean(irToAverage(:, :, validReps), 3);
        else
            irAvg = zeros(Lh, numErrMics);
            fprintf('[ERROR] Spk%d 所有重复均被丢弃！\n', spk);
        end
        impulseResponses(:, :, spk) = irAvg;

        if cfg.exportAlignedIR
            impulseResponsesAligned(:, :, spk) = irAvg;
        end

        % === 6. 延迟推荐 & 可用性判定 ===
        driftStable = ~isempty(allReliablePeaks) && ...
            (iqr(allReliablePeaks) <= cfg.reliableMaxIQR) && ...
            (max(allReliablePeaks) - min(allReliablePeaks) <= cfg.maxAllowedDriftSamples);

        medianSNR_perMic = zeros(1, numErrMics);
        for m = 1:numErrMics
            micSNR = repSNRlist(validReps, m);
            validSNR = micSNR(isfinite(micSNR) & ~isnan(micSNR));
            if ~isempty(validSNR)
                medianSNR_perMic(m) = median(validSNR);
            else
                medianSNR_perMic(m) = -Inf;
            end
        end
        medianSNR = median(medianSNR_perMic);
        
        % 修复：处理 median 计算中的无效语法
        medianEnergyConcentration = NaN;
        medianPeakSharpness = NaN;
        if any(validReps)
            tempEC = repEnergyConcentration(validReps, :);
            medianEnergyConcentration = median(tempEC(:), 'omitnan');
            
            tempPS = repPeakSharpness(validReps, :);
            medianPeakSharpness = median(tempPS(:), 'omitnan');
        end

        reliableRatioAll = 0;
        if any(validReps)
            % 修复：正确连接所有有效重复的可靠掩码
            validRelCells = repPeakReliable(validReps);
            if ~isempty(validRelCells)
                validRelMask = [];
                for i = 1:length(validRelCells)
                    validRelMask = [validRelMask; validRelCells{i}(:)];
                end
                if ~isempty(validRelMask)
                    numValidReliable = sum(validRelMask);
                    numValidTotal = sum(validReps) * numErrMics;
                    reliableRatioAll = numValidReliable / numValidTotal;
                end
            end
        end

        if isempty(allReliablePeaks)
            [~, peakPos] = max(abs(irAvg), [], 1);
            recommendedDelayReliable = round(median(peakPos));
            fprintf('[WARN] Spk%d 没有可靠峰值，使用幅度最大位置: %d\n', spk, recommendedDelayReliable);
        else
            recommendedDelayReliable = round(median(allReliablePeaks));
        end

        usable = driftStable && ...
             (reliableRatioAll >= cfg.reliableMinRatio) && ...
             (~isempty(allReliablePeaks)) && ...
             (medianSNR >= cfg.snrThresholdDB);

        meta.perSpeaker{spk} = struct(...
            'repPeakPositions', repPeakLists_raw, ...
            'repPeakReliable', repPeakReliable, ...
            'repSNRlist', repSNRlist, ...
            'repEnergyConcentration', repEnergyConcentration, ...
            'repPeakSharpness', repPeakSharpness, ...
            'driftStable', driftStable, ...
            'reliableRatioAll', reliableRatioAll, ...
            'medianSNR', medianSNR, ...
            'medianEnergyConcentration', medianEnergyConcentration, ...
            'medianPeakSharpness', medianPeakSharpness, ...
            'recommendedDelayReliable', recommendedDelayReliable, ...
            'usable', usable ...
            );

        fprintf('[diag] Spk%d 可用性=%d, medianSNR=%.2f dB, 可靠比例=%.1f%%, Energy85=%.0f, PeakSharpness=%.1f\n', ...
            spk, usable, medianSNR, reliableRatioAll*100, medianEnergyConcentration, medianPeakSharpness);
        
        if medianSNR < cfg.snrThresholdDB && cfg.enableAutoGainAdjustment && spk < cfg.numSpeakers
            fprintf('[diag] 警告: Spk%d SNR过低 (%.1f dB < %.1f dB)，建议:\n', ...
                spk, medianSNR, cfg.snrThresholdDB);
            fprintf('       1. 提高spkAmplitude (当前: %.2f)\n', cfg.spkAmplitude(spk));
            fprintf('       2. 检查麦克风位置和连接\n');
            fprintf('       3. 降低环境噪声\n');
        end
    end
        
    %% ============== 构建输出结构 ==============
    secondary = struct();
    secondary.impulseResponses = impulseResponses;
    if cfg.exportAlignedIR
        secondary.impulseResponses_aligned = impulseResponsesAligned;
    end
    secondary.fs = cfg.fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;
    secondary.errorMicPhysicalChannels = errMicIdx;
    secondary.irLength = Lh;
    secondary.description = 'Secondary path v3.2 (fixed alignment)';
    secondary.timestampUtc = datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS');
    secondary.measurementConfig = cfg;

    recommendedVector = arrayfun(@(i) meta.perSpeaker{i}.recommendedDelayReliable, 1:cfg.numSpeakers);
    secondary.delayEstimateSamples = recommendedVector;
    usableVector = arrayfun(@(i) meta.perSpeaker{i}.usable, 1:cfg.numSpeakers);
    secondary.usable = usableVector;
    
    if cfg.saveDiagnosticInfo
        secondary.meta = meta;
        if cfg.saveEachRepeatIR
            secondary.irEachRepeat = irEachRepeat;
        end
        if saveRaw
            secondary.rawRecordings = rawRecordings;
        end
    end

    mkdir(fileparts(cfg.secondaryPathFile));
    save(cfg.secondaryPathFile, 'secondary', '-v7.3');
    fprintf('[measure] 保存完成 -> %s\n', cfg.secondaryPathFile);
    
    generate_measurement_report(secondary, meta);

catch ME
    fprintf('[ERROR] %s\n', ME.message);
    fprintf('[ERROR] Stack trace:\n');
    for k = 1:length(ME.stack)
        fprintf('  %s:%d\n', ME.stack(k).name, ME.stack(k).line);
    end

    % 安全释放硬件（仅当有效时）
    if exist('hw', 'var') && ~isempty(hw) && ismethod(hw, 'release')
        try
            hw.release();
        catch relErr
            fprintf('[WARN] hw.release() failed: %s\n', relErr.message);
        end
    end
    rethrow(ME); % 重新抛出异常以便调试
end

fprintf('[measure] 测量完成，总计 %d 个扬声器，%d 个误差麦克风。\n', ...
    cfg.numSpeakers, numErrMics);