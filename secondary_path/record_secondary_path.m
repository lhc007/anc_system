% record_secondary_path.m 多扬声器→多误差麦克风次级路径测量（ESS）
% 优化点：多麦对齐、IR峰值对齐、内存控制、动态SNR、错误安全释放

clear; clc;

%% ============== 配置区 ==============
cfg = anc_config();

% 仅处理误差麦克风（关键！避免处理无用通道）
errMicIdx = cfg.micChannels.error;               % e.g., [5 6]
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
    'fftCorrEnable', cfg.deconvFftCorrEnable ...
    );

%% ============== 初始化硬件（带安全释放）=============
try


    %% ============== 生成 Sweep ==============
    fprintf('[measure] 生成 ESS sweep...\n');
    [sweepSigCore, invFilter, sweepInfo] = generate_sweep(sweepCfg);
    preSilence  = zeros(round(cfg.preSilenceSec*cfg.fs),1);
    postSilence = zeros(round(cfg.postSilenceSec*cfg.fs),1);

    if cfg.enableLowFreqBoost
        [b_lp,a_lp] = butter(4, cfg.lowFreqCutHz/(cfg.fs/2), 'low');
        sweepLow    = filter(b_lp,a_lp,sweepSigCore);
        sweepMixed  = (1 - cfg.lowFreqMixRatio)*sweepSigCore + cfg.lowFreqMixRatio*sweepLow;
        sweepMixed  = sweepMixed / max(abs(sweepMixed)+1e-12) * sweepCfg.amplitude;
        sweepDrive  = sweepMixed;
    else
        sweepDrive  = sweepSigCore / max(abs(sweepSigCore)+1e-12) * sweepCfg.amplitude;
    end

    sweepSig = [preSilence; sweepDrive; postSilence];
    totalSamples = length(sweepSig);

    % ✅ 新增：记录核心 sweep 的起止位置
    sweepCoreStart = length(preSilence) + 1;
    sweepCoreEnd   = sweepCoreStart + length(sweepDrive) - 1;
    sweepCoreLen   = length(sweepDrive);
    fprintf('[measure] Sweep 总长=%d 样本\n', totalSamples);

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
    meta = struct();
    meta.global = struct('sampleRate',cfg.fs,'sweepLen',totalSamples,'repetitions',cfg.repetitions,...
        'minPhysDelaySamples',cfg.minPhysDelaySamples,'enableRepeatAlignment',cfg.enableRepeatAlignment,...
        'lowFreqBoost',cfg.enableLowFreqBoost,'exportAlignedIR',cfg.exportAlignedIR);
    meta.perSpeaker = cell(cfg.numSpeakers,1);
    allIrRepList = cell(cfg.numSpeakers,1);
    %% ============== 测量循环 ==============
    hw = [];
    for spk = 1:cfg.numSpeakers
        fprintf('\n[measure] ===== 扬声器 %d =====\n', spk);

        repGlobalShifts = zeros(cfg.repetitions,1);
        repRecordedRaw = cell(cfg.repetitions,1);
        repDiscardRepeat = false(cfg.repetitions,1);
        repPeakLists_raw = cell(cfg.repetitions,1);
        repPeakLists_forIR = cell(cfg.repetitions,1);
        repPeakReliable = cell(cfg.repetitions,1);
        repSNRlist = zeros(cfg.repetitions, numErrMics);
        irRepList = zeros(Lh, numErrMics, cfg.repetitions);

        % === 阶段1：采集所有 repeats ===
        for rep = 1:cfg.repetitions
            fprintf('[measure]  播放重复 %d/%d...\n', rep, cfg.repetitions);
            
            if rep == 1
            fprintf('[measure] 初始化硬件...\n');
            hw = hardware_init_measure(cfg);
            else
                % 释放上一次的硬件
                hw.release();
                pause(0.3);  % 给 Windows 音频子系统时间清理资源
                hw = hardware_init_measure(cfg);  % 重新打开
            end

            % === 连续指针模式播放 ===
            spkDriveSig = cfg.spkAmplitude(spk) * sweepSig(:);  % 列向量

            % 防削波（仅针对该扬声器信号）
            maxVal = max(abs(spkDriveSig));
            if maxVal > 0.98
                scale = 0.95 / maxVal;
                spkDriveSig = spkDriveSig * scale;
                fprintf('[measure] 播放缩放 %.2f\n', scale);
            end
            % 保存用于 xcorr 的参考信号（含实际缩放）
            refSweepForXcorr = spkDriveSig;  % 长度 = totalSamples

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
                % 构造输出块（全零初始化）
                outBlock = zeros(blockSize, cfg.numSpeakers);
                % 从 sweep 连续读取
                if ptr <= totalSamples
                    nToCopy = min(blockSize, totalSamples - ptr + 1);
                    outBlock(1:nToCopy, spk) = spkDriveSig(ptr : ptr + nToCopy - 1);
                    ptr = ptr + nToCopy;
                end
                % 播放并录制
                hw.writer(outBlock);
                micFrame = hw.reader();

                % 安全处理 micFrame 尺寸
                if isempty(micFrame)
                    micFrame = zeros(blockSize, cfg.micNumChannels);
                elseif size(micFrame, 1) < blockSize
                    micFrame(end+1:blockSize, :) = 0;
                end
                recordedFull(rdPtr : rdPtr + blockSize - 1, :) = micFrame;
                rdPtr = rdPtr + blockSize;
            end

            % 截取与原始 sweep 等长（确保后续 xcorr 对齐）
            recordedFull = recordedFull(1:totalSamples, :);
            recorded = recordedFull(:, errMicIdx);

            if saveRaw && ((spk==1 && rep==1 && cfg.saveFirstRaw) || cfg.saveAllRaw)
                rawRecordings{spk,rep} = recorded;
            end
            repRecordedRaw{rep} = recorded;

            % === 改进：使用多误差麦加权互相关估计 shift ===
            shiftEst = 0;
            weightSum = 0;
            for m = 1:numErrMics
                [c, lags] = xcorr(recorded(:,m), refSweepForXcorr, 'coeff');

                % [~, idx] = max(c); lag = lags(idx);
                validRange = (lags >= -cfg.fs*0.5) & (lags <= cfg.fs*2);  % 示例：-0.5s ~ +2s
                if any(validRange)
                    [~, idx] = max(c(validRange));
                    lag = lags(validRange);
                    lag = lag(idx);
                else
                    lag = NaN;
                end
                w = max(abs(c));  % 相关强度作权重
                shiftEst = shiftEst + w * lag;
                weightSum = weightSum + w;
            end
            globalShift = round(shiftEst / (weightSum + 1e-12));
            repGlobalShifts(rep) = globalShift;

            fprintf('[measure]  Spk%d Rep%d globalShift=%d\n', spk, rep, globalShift);
        end
        % 所有 repeats 结束后释放
        hw.release();
        % === 阶段1.5：确定 refDelay（✅ 修复：使用 globalShift 中位数作为真实系统延迟）===
        validShifts = repGlobalShifts(~isnan(repGlobalShifts) & isfinite(repGlobalShifts));
        if ~isempty(validShifts)
            refDelay = round(median(validShifts));
            % ✅ 新增：强制 refDelay 不小于物理最小延迟
            refDelay = max(refDelay, cfg.minPhysDelaySamples);
        else
            refDelay = cfg.minPhysDelaySamples;  % fallback to physical min
        end
        fprintf('[measure]  Spk%d refDelay set to %d samples (median of globalShift)\n', spk, refDelay);

        alignThreshold = max(cfg.minPhysDelaySamples/2, cfg.maxAllowedDriftSamples);
        for rep = 1:cfg.repetitions
            if isnan(repGlobalShifts(rep)) || abs(repGlobalShifts(rep)-refDelay) > alignThreshold
                repDiscardRepeat(rep) = true;
                fprintf('[warn]  Spk%d Rep%d 剔除 (shift=%d)\n', spk, rep, repGlobalShifts(rep));
            end
        end

        % === 阶段2：对齐 + 反卷积 ===
        for rep = 1:cfg.repetitions
            rec = repRecordedRaw{rep};
            shiftSamples = repGlobalShifts(rep) - refDelay;
            if shiftSamples > 0
                recShift = [rec(shiftSamples+1:end,:); zeros(shiftSamples, numErrMics)];
            elseif shiftSamples < 0
                s = -shiftSamples;
                recShift = [zeros(s, numErrMics); rec(1:end-s,:)];
            else
                recShift = rec;
            end

            irCurrent = zeros(Lh, numErrMics);
            peakIdxEach_raw = zeros(numErrMics, 1);
            peakRelEach = false(numErrMics, 1);
            snrMic = zeros(numErrMics, 1);

            for m = 1:numErrMics
                % 至少覆盖 sweep + 最大 IR 长度
                L_needed = sweepCoreLen + cfg.irMaxLen;
                L_extract = min(L_needed, size(recShift,1));
                if L_extract < sweepCoreLen + cfg.minPhysDelaySamples
                    warning('[measure] Spk%d Rep%d: recShift too short for reliable deconv (%d < %d)', ...
                        spk, rep, L_extract, sweepCoreLen + cfg.minPhysDelaySamples);
                end
                recForDeconv = recShift(1:L_extract, m);
                % 不含静音的核心 sweep
                sweepForDeconv = sweepDrive;

                outStruct = deconvolve_sweep(recForDeconv, sweepForDeconv, cfg.fs, deconvParams);
                h_full = outStruct.h;
                pk = outStruct.peakIdx;
                snrMic(m) = outStruct.snrEst;

                % 改进：使用 cumulative energy 定位主峰（更鲁棒）
                if pk < cfg.minPhysDelaySamples
                    % cumE = cumsum(abs(h_full).^2);
                    energySeq = abs(h_full).^2;
                    totalE = sum(energySeq);

                    if totalE < 1e-18  % 判定为全零或近零
                        newPk = [];    % 无有效峰值
                    else
                        cumE = cumsum(energySeq) / totalE;      % 归一化，cumE(end) == 1
                        targetE = cfg.deconvCumEnergyFrac;      % e.g., 0.95
                        newPk = find(cumE >= targetE, 1, 'first');
                    end
                    if ~isempty(newPk) && newPk > pk
                        pk = newPk;
                    end
                end

                reliable = (pk >= cfg.minPhysDelaySamples) && ...
                    (outStruct.snrEst >= cfg.snrThresholdDB) && ...
                    (isfield(outStruct, 'peakReliability') && outStruct.peakReliability);

                peakIdxEach_raw(m) = pk;
                peakRelEach(m) = reliable;
                Luse = min(Lh, length(h_full));
                irCurrent(1:Luse, m) = h_full(1:Luse);
            end

            repPeakLists_raw{rep} = peakIdxEach_raw;
            repPeakReliable{rep} = peakRelEach;
            repSNRlist(rep, :) = snrMic;
            irRepList(:, :, rep) = irCurrent;
        end

        % === 可用性判定（仅误差麦）===
        validRepMask = ~repDiscardRepeat;
        reliablePeaksAll = [];
        for rep = 1:cfg.repetitions
            if ~validRepMask(rep), continue; end
            reliablePeaksAll = [reliablePeaksAll; repPeakLists_raw{rep}(repPeakReliable{rep})];
        end

        driftStable = ~isempty(reliablePeaksAll) && ...
            (iqr(reliablePeaksAll) <= cfg.reliableMaxIQR) && ...
            (max(reliablePeaksAll)-min(reliablePeaksAll) <= cfg.maxAllowedDriftSamples);

        % 可靠率 & median SNR
        if any(validRepMask)
            numValidReliable = sum([repPeakReliable{validRepMask}]);
            numValidTotal = sum(validRepMask) * numErrMics;
            reliableRatioAll = numValidReliable / max(1, numValidTotal);
            medianSNR = median(repSNRlist(validRepMask, :), 'all');
        else
            reliableRatioAll = 0; medianSNR = -Inf;
        end

        % 平均 IR
        irAvgStruct = average_ir(irRepList, cfg.doAlignRepeats);
        irAvg = irAvgStruct.irAvg;
        impulseResponses(:, :, spk) = irAvg;

        if cfg.exportAlignedIR
            irAlignedStruct = average_ir(irRepList, true);
            impulseResponsesAligned(:, :, spk) = irAlignedStruct.irAvg;
        end

        % 推荐延迟
        if isempty(reliablePeaksAll)
            % 若无可靠峰值，用 IR 幅度最大值位置的中位数作为 fallback
            [~, peakPos] = max(abs(irAvg), [], 1);
            recommendedDelayReliable = round(median(peakPos));
        else
            recommendedDelayReliable = round(median(reliablePeaksAll));
        end

        % 可用性判定（含 SNR 门限）
        usable = driftStable && ...
            (reliableRatioAll >= cfg.reliableMinRatio) && ...
            (~isempty(reliablePeaksAll)) && ...
            (medianSNR >= cfg.snrThresholdDB);

        % 保存元数据
        meta.perSpeaker{spk} = struct(...
            'repGlobalShifts', repGlobalShifts, ...
            'driftStable', driftStable, ...
            'reliableRatioAll', reliableRatioAll, ...
            'medianSNR', medianSNR, ...
            'recommendedDelayReliable', recommendedDelayReliable, ...
            'usable', usable, ...
            'repDiscardRepeat', repDiscardRepeat ...
            );

        fprintf('[diag] Spk%d usable=%d medianSNR=%.2f dB\n', spk, usable, medianSNR);
    end

    %% ============== 构建输出结构 ==============
    secondary = struct();
    secondary.impulseResponses = impulseResponses;  % [L, numErrMics, numSpeakers]
    if cfg.exportAlignedIR
        secondary.impulseResponses_aligned = impulseResponsesAligned;
    end
    secondary.fs = cfg.fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;  % ⚠️ 关键：只记录误差麦数量
    secondary.errorMicPhysicalChannels = errMicIdx;  % 明确物理通道
    secondary.irLength = Lh;
    secondary.description = 'Secondary path v2.5 (optimized)';
    secondary.timestampUtc = datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS');

    % 推荐延迟向量
    recommendedVector = arrayfun(@(i) meta.perSpeaker{i}.recommendedDelayReliable, 1:cfg.numSpeakers);
    secondary.delayEstimateSamples = recommendedVector;

    % 保存
    mkdir(fileparts(cfg.secondaryPathFile));
    save(cfg.secondaryPathFile, 'secondary', '-v7.3');
    fprintf('[measure] 保存完成 -> %s\n', cfg.secondaryPathFile);

catch ME
    fprintf('[ERROR] %s\n', ME.message);
    if ~isempty(hw)
        try
            hw.release();
        catch
            fprintf('[ERROR] 释放资源错误 %s\n', ME.message);
        end
    end
end

fprintf('[measure] Done.\n');