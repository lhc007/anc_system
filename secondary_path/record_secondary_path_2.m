% record_secondary_path.m 多扬声器→多误差麦克风次级路径测量（ESS）
% 优化点：多麦对齐、IR峰值对齐、内存控制、动态SNR、错误安全释放
% 版本：v2.6 优化版

clear; clc;

%% ============== 配置区 ==============
cfg = anc_config();

% 仅处理误差麦克风（关键！避免处理无用通道）
errMicIdx = cfg.micChannels.error;               % e.g., [5 6]
numErrMics = length(errMicIdx);
fprintf('[measure] 误差麦克风通道: %s\n', mat2str(errMicIdx));

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
    [sweepSigCore, invFilter, sweepInfo] = generate_sweep(cfg);
    preSilence  = zeros(round(cfg.preSilenceSec*cfg.fs),1);
    postSilence = zeros(round(cfg.postSilenceSec*cfg.fs),1);

    if cfg.enableLowFreqBoost
        [b_lp,a_lp] = butter(4, cfg.lowFreqCutHz/(cfg.fs/2), 'low');
        sweepLow    = filter(b_lp,a_lp,sweepSigCore);
        sweepMixed  = (1 - cfg.lowFreqMixRatio)*sweepSigCore + cfg.lowFreqMixRatio*sweepLow;
        sweepMixed  = sweepMixed / max(abs(sweepMixed)+1e-12) * cfg.amplitude;
        sweepDrive  = sweepMixed;
    else
        sweepDrive  = sweepSigCore / max(abs(sweepSigCore)+1e-12) * cfg.amplitude;
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
    
    % 新增：保存每个repeat的IR用于诊断
    if cfg.saveEachRepeatIR
        irEachRepeat = cell(cfg.numSpeakers, 1);
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
        repEnergyConcentration = zeros(cfg.repetitions, numErrMics); % 新增：能量集中度
        repPeakSharpness = zeros(cfg.repetitions, numErrMics); % 新增：峰值锐度
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

        % === 阶段2：基于 globalShift 窗口截取 + 反卷积（不依赖 shift 对齐）===
        % 确定参考延迟窗口起点（仅用于截取）- 使用trimmed mean提高鲁棒性
        validShifts = repGlobalShifts(~isnan(repGlobalShifts) & isfinite(repGlobalShifts));
        if ~isempty(validShifts)
            if length(validShifts) > 2
                % 使用trimmed mean：去掉最小和最大值
                validShiftsSorted = sort(validShifts);
                trimmed = validShiftsSorted(2:end-1);  % 去掉最小和最大
                refWindowCenter = round(mean(trimmed));
            else
                refWindowCenter = round(median(validShifts));
            end
        else
            refWindowCenter = cfg.minPhysDelaySamples + sweepCoreLen/2; % fallback
        end

        % 定义安全窗口：确保覆盖最小物理延迟到最大可能 IR 结束
        winStartMin = max(1, cfg.minPhysDelaySamples);
        winEndMax   = sweepCoreLen + cfg.irMaxLen;
        winHalfLen  = ceil(cfg.irMaxLen / 2);  % 窗口半长，围绕估计峰值

        for rep = 1:cfg.repetitions
            rec = repRecordedRaw{rep};
            globalShift = repGlobalShifts(rep);

            % === 1. 用 globalShift 定位窗口中心（仅用于截取）===
            if isnan(globalShift) || ~isfinite(globalShift)
                winCenter = refWindowCenter;
            else
                winCenter = globalShift;
            end
            winCenter = max(winStartMin + winHalfLen, min(winEndMax - winHalfLen, winCenter));

            % 截取窗口：[winCenter - winHalfLen, winCenter + winHalfLen + sweepCoreLen]
            extractStart = max(1, winCenter - winHalfLen);
            extractEnd   = min(size(rec,1), extractStart + sweepCoreLen + cfg.irMaxLen - 1);
            
            % ✅ 新增：检查窗口长度是否足够
            expectedLen = sweepCoreLen + cfg.irMaxLen;
            actualLen = extractEnd - extractStart + 1;
            if actualLen < expectedLen
                warning('[measure] Spk%d Rep%d: 窗口长度不足 (实际:%d < 预期:%d)', ...
                    spk, rep, actualLen, expectedLen);
            end
            
            recExtract   = rec(extractStart:extractEnd, :);
            offsetInRec  = extractStart - 1;  % 用于将 IR 峰值映射回全局时间轴

            % === 2. 对每个误差麦反卷积 ===
            irCurrent = zeros(Lh, numErrMics);
            peakIdxEach_raw = zeros(numErrMics, 1);
            peakRelEach = false(numErrMics, 1);
            snrMic = zeros(numErrMics, 1);
            energyConcentration = zeros(numErrMics, 1);
            peakSharpness = zeros(numErrMics, 1);

            for m = 1:numErrMics
                outStruct = deconvolve_sweep(recExtract(:,m), sweepDrive, cfg.fs, deconvParams);
                h_full = outStruct.h;
                pk_local = outStruct.peakIdx;  % 峰值在 h_full 中的位置

                % 映射回全局样本索引（相对于原始 recordedFull 起点）
                pk_global = offsetInRec + pk_local;

                % ✅ 新增：IR质量指标计算
                if ~isempty(h_full)
                    % 1. 能量集中度（85%能量位置）
                    energySeq = abs(h_full).^2;
                    totalE = sum(energySeq);
                    if totalE > 1e-18
                        cumE = cumsum(energySeq) / totalE;
                        energy85 = find(cumE >= 0.85, 1, 'first');
                        energyConcentration(m) = energy85;
                        
                        % 2. 峰值锐度（峰值/有效RMS）
                        ir_abs = abs(h_full);
                        peakVal = max(ir_abs);
                        % 计算有效区域的RMS（排除前5个样本的预延迟）
                        validStart = max(1, min(length(h_full), cfg.minPhysDelaySamples+1));
                        validEnd = min(length(h_full), validStart + cfg.irMaxLen - 1);
                        if validEnd >= validStart
                            validIR = h_full(validStart:validEnd);
                            rmsValid = sqrt(mean(abs(validIR).^2));
                            if rmsValid > 0
                                peakSharpness(m) = peakVal / rmsValid;
                            end
                        end
                    end
                end

                % 改进：若主峰太早，尝试用 cumulative energy 修正
                if pk_local < cfg.minPhysDelaySamples && ~isempty(h_full)
                    energySeq = abs(h_full).^2;
                    totalE = sum(energySeq);
                    if totalE > 1e-18
                        cumE = cumsum(energySeq) / totalE;
                        newPk_local = find(cumE >= cfg.deconvCumEnergyFrac, 1, 'first');
                        if ~isempty(newPk_local) && newPk_local > pk_local
                            pk_local = newPk_local;
                            pk_global = offsetInRec + pk_local;
                        end
                    end
                end

                 % ========== 【新增】DEBUG 打印 ==========
                snrEst_val = outStruct.snrEst;
                hasPeakRel = isfield(outStruct, 'peakReliability');
                peakRel_val = hasPeakRel && outStruct.peakReliability;
                
                delay_ok = (pk_global >= cfg.minPhysDelaySamples);
                snr_ok = (snrEst_val >= cfg.snrThresholdDB);
                
                delay_status = 'FAIL'; if delay_ok, delay_status = 'OK'; end
                snr_status   = 'FAIL'; if snr_ok,   snr_status   = 'OK'; end
                
                reliable_val = delay_ok && snr_ok && peakRel_val;
                
                fprintf('  [DEBUG] Mic%d: pk_global=%d (minDelay=%d → %s), SNR=%.2f dB (thresh=%.1f → %s), peakRel=%d → reliable=%d\n', ...
                    m, pk_global, cfg.minPhysDelaySamples, delay_status, ...
                    snrEst_val, cfg.snrThresholdDB, snr_status, ...
                    peakRel_val, reliable_val);
                % =======================================

                reliable = (pk_global >= cfg.minPhysDelaySamples) && ...
                           (outStruct.snrEst >= cfg.snrThresholdDB) && ...
                           (isfield(outStruct, 'peakReliability') && ...
                           outStruct.peakReliability);
                
                peakIdxEach_raw(m) = pk_global;
                peakRelEach(m) = reliable;
                snrMic(m) = outStruct.snrEst;
                
                % 保存质量指标
                energyConcentration(m) = energyConcentration(m);
                peakSharpness(m) = peakSharpness(m);

                Luse = min(Lh, length(h_full));
                irCurrent(1:Luse, m) = h_full(1:Luse);
            end

            repPeakLists_raw{rep} = peakIdxEach_raw;
            repPeakReliable{rep} = peakRelEach;
            repSNRlist(rep, :) = snrMic;
            repEnergyConcentration(rep, :) = energyConcentration;
            repPeakSharpness(rep, :) = peakSharpness;
            irRepList(:, :, rep) = irCurrent;
        end
        
        % 保存每个repeat的IR用于后期分析
        if cfg.saveEachRepeatIR
            irEachRepeat{spk} = irRepList;
        end

        % === 3. 基于 IR 峰值做重复间对齐（若启用）===
        % 先收集所有可靠峰值用于对齐参考
        allReliablePeaks = [];
        for rep = 1:cfg.repetitions
            peaks = repPeakLists_raw{rep}(repPeakReliable{rep});
            allReliablePeaks = [allReliablePeaks; peaks];
        end

        if isempty(allReliablePeaks)
            refPeakForAlign = cfg.minPhysDelaySamples + 50; % fallback
        else
            refPeakForAlign = round(median(allReliablePeaks));
        end

        % 如果需要对齐 repeats（用于平均）
        if cfg.doAlignRepeats || cfg.exportAlignedIR
            irRepList_aligned = zeros(Lh, numErrMics, cfg.repetitions);
            for rep = 1:cfg.repetitions
                ir_orig = irRepList(:, :, rep);
                % 计算该 repeat 的平均可靠峰值（若无可靠，则用整体 median）
                reliableMask = repPeakReliable{rep};
                if any(reliableMask)
                    repPeak = mean(repPeakLists_raw{rep}(reliableMask));
                else
                    repPeak = refPeakForAlign;
                end
                shiftIR = round(refPeakForAlign - repPeak);  % 正值表示需右移（补零在前）

                for m = 1:numErrMics
                    h = ir_orig(:, m);
                    if shiftIR > 0
                        h_shifted = [zeros(shiftIR,1); h(1:end-shiftIR)];
                    elseif shiftIR < 0
                        h_shifted = [h(-shiftIR+1:end); zeros(-shiftIR,1)];
                    else
                        h_shifted = h;
                    end
                    irRepList_aligned(:, m, rep) = h_shifted(1:Lh);
                end
            end
        else
            irRepList_aligned = irRepList;
        end

        % === 4. 平均 IR ===
        irAvg = mean(irRepList_aligned, 3);
        impulseResponses(:, :, spk) = irAvg;

        if cfg.exportAlignedIR
            impulseResponsesAligned(:, :, spk) = irAvg; % 已对齐
        end

        % === 5. 延迟推荐 & 可用性判定（完全基于 IR 峰值）===
        reliablePeaksAll = allReliablePeaks; % 已收集

        driftStable = ~isempty(reliablePeaksAll) && ...
            (iqr(reliablePeaksAll) <= cfg.reliableMaxIQR) && ...
            (max(reliablePeaksAll) - min(reliablePeaksAll) <= cfg.maxAllowedDriftSamples);

        % ✅【修复】优化medianSNR计算：先按麦克风计算中位数，再整体中位数
        medianSNR_perMic = zeros(1, numErrMics);
        for m = 1:numErrMics
            micSNR = repSNRlist(:, m);
            validSNR = micSNR(isfinite(micSNR) & ~isnan(micSNR));
            if ~isempty(validSNR)
                medianSNR_perMic(m) = median(validSNR);
            else
                medianSNR_perMic(m) = -Inf;
            end
        end
        medianSNR = median(medianSNR_perMic);

        % 计算IR质量指标的中位数
        medianEnergyConcentration = median(repEnergyConcentration(:));
        medianPeakSharpness = median(repPeakSharpness(:));

        % 可靠性比例（用于 usable 判定）
        validRelMask = [repPeakReliable{:}];
        if any(validRelMask)
            numValidReliable = sum(validRelMask);
            numValidTotal = cfg.repetitions * numErrMics;
            reliableRatioAll = numValidReliable / numValidTotal;
        else
            reliableRatioAll = 0;
        end

        if isempty(reliablePeaksAll)
            % fallback：用平均 IR 的幅度最大值
            [~, peakPos] = max(abs(irAvg), [], 1);
            recommendedDelayReliable = round(median(peakPos));
        else
            recommendedDelayReliable = round(median(reliablePeaksAll));
        end

        % 最终可用性判定（严格）
        usable = driftStable && ...
                 (reliableRatioAll >= cfg.reliableMinRatio) && ...
                 (~isempty(reliablePeaksAll)) && ...
                 (medianSNR >= cfg.snrThresholdDB);

        % 保存元数据（包含新增的质量指标）
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

        fprintf('[diag] Spk%d usable=%d medianSNR=%.2f dB, Energy85=%.0f, PeakSharpness=%.1f\n', ...
            spk, usable, medianSNR, medianEnergyConcentration, medianPeakSharpness);
        
        % ✅ 新增：如果SNR太低但可重测，建议调整增益
        if medianSNR < cfg.snrThresholdDB && cfg.enableAutoGainAdjustment && spk < cfg.numSpeakers
            fprintf('[diag] 警告: Spk%d SNR过低 (%.1f dB < %.1f dB)，考虑提高spkAmplitude或检查连接\n', ...
                spk, medianSNR, cfg.snrThresholdDB);
        end
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
    secondary.description = 'Secondary path v2.6 (optimized with quality metrics)';
    secondary.timestampUtc = datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS');

    % 推荐延迟向量
    recommendedVector = arrayfun(@(i) meta.perSpeaker{i}.recommendedDelayReliable, 1:cfg.numSpeakers);
    secondary.delayEstimateSamples = recommendedVector;
    
    % 保存诊断信息
    if cfg.saveDiagnosticInfo
        secondary.meta = meta;
        if cfg.saveEachRepeatIR
            secondary.irEachRepeat = irEachRepeat;
        end
    end

    % 保存
    mkdir(fileparts(cfg.secondaryPathFile));
    save(cfg.secondaryPathFile, 'secondary', '-v7.3');
    fprintf('[measure] 保存完成 -> %s\n', cfg.secondaryPathFile);

catch ME
    fprintf('[ERROR] %s\n', ME.message);
    fprintf('[ERROR] Stack trace:\n');
    for k = 1:length(ME.stack)
        fprintf('  %s:%d\n', ME.stack(k).name, ME.stack(k).line);
    end
    
    if ~isempty(hw)
        try
            hw.release();
        catch
            fprintf('[ERROR] 释放资源错误\n');
        end
    end
    rethrow(ME);
end

fprintf('[measure] Done.\n');