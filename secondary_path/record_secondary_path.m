% record_secondary_path_v2_4.m
% 多扬声器→多麦克风次级路径测量 (ESS)
% - 支持 ASIO4ALL 聚合多声卡为单设备（4通道输出）
% - 漂移补偿 / 可靠峰统计 / 双延迟推荐 / 可用性判定 / 对齐IR
% - v2.4: 重构为单 writer 架构，解决跨设备异步问题

clear; clc;

%% ============== 配置区 ==============
cfg = anc_config();

% Sweep
sweepCfg = struct('fs',cfg.fs,'T',2.0,'f1',20,'f2',6000,...
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
    'envSmoothWin',  cfg.envSmoothWin , ...
    'cumEnergyFrac', cfg.deconvCumEnergyFrac, ...
    'minPeakFrac',   cfg.deconvMinPeakFrac, ...
    'snrBodyRadius', cfg.deconvSnrBodyRadius, ...
    'fftCorrEnable', cfg.deconvFftCorrEnable ...
);

fprintf('[measure] 初始化硬件...\n');
hw = hardware_init_measure(cfg);  % ← 现在只返回 reader + writer（4ch）

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
fprintf('[measure] Sweep 总长=%d 样本\n', length(sweepSig));

%% ============== 容器 ==============
Lh = cfg.irMaxLen;
impulseResponses        = zeros(Lh,cfg.micNumChannels,cfg.numSpeakers);
impulseResponsesAligned = zeros(Lh,cfg.micNumChannels,cfg.numSpeakers);

blockSize    = cfg.timeFrameSamples;
totalSamples = length(sweepSig);
numBlocks    = ceil(totalSamples/blockSize);

rawRecordings = cell(cfg.numSpeakers,cfg.repetitions);

meta = struct();
meta.global = struct('sampleRate',cfg.fs,'sweepLen',totalSamples,'repetitions',cfg.repetitions,...
    'minPhysDelaySamples',cfg.minPhysDelaySamples,'enableRepeatAlignment',cfg.enableRepeatAlignment,...
    'lowFreqBoost',cfg.enableLowFreqBoost,'exportAlignedIR',cfg.exportAlignedIR);
meta.perSpeaker = cell(cfg.numSpeakers,1);
allIrRepList = cell(cfg.numSpeakers,1);

%% ============== 测量循环 ==============
for spk=1:cfg.numSpeakers
    fprintf('\n[measure] ===== 扬声器 %d =====\n', spk);
    irRepList = zeros(Lh,cfg.micNumChannels,cfg.repetitions);

    repDelaysRaw      = zeros(cfg.repetitions,1);
    repDelaysCorr     = zeros(cfg.repetitions,1);
    repPeakLists      = cell(cfg.repetitions,1);
    repPeakReliable   = cell(cfg.repetitions,1);
    repSNRlist        = zeros(cfg.repetitions,cfg.micNumChannels);
    repDiscardRepeat  = false(cfg.repetitions,1);

    refDelay = NaN;

    for rep=1:cfg.repetitions
        fprintf('[measure]  播放重复 %d/%d...\n', rep,cfg.repetitions);
        
        % === 构建 4 通道激励信号（仅当前扬声器有信号）===
        outAll = zeros(totalSamples, 4);
        ampVec = cfg.spkAmplitude(:)'; % [1.0, 1.0, 1.5, 1.5]
        for ch = 1:4
            if spk == ch
                outAll(:, ch) = ampVec(ch) * sweepSig;
            end
        end
        
        % 归一化防削波
        maxPlay = max(abs(outAll(:)));
        if maxPlay > 0.98
            scale = 0.95 / maxPlay;
            outAll = outAll * scale;
            fprintf('[measure] 播放缩放 %.2f\n', scale);
        end

        recorded = zeros(numBlocks*blockSize, cfg.micNumChannels);
        
        % 预热：发送静音
        for pr = 1:cfg.preRollFrames
            hw.writer(zeros(blockSize, 4));  % ← 4 通道零信号
            hw.reader();
        end
        
        wrPtr = 1;
        for b = 1:numBlocks
            iStart = (b-1)*blockSize + 1;
            iEnd = min(b*blockSize, totalSamples);
            blk = outAll(iStart:iEnd, :);
            if cfg.writeBlockPad && size(blk,1) < blockSize
                blk(end+1:blockSize, :) = 0;
            end
            hw.writer(blk);  % ← 单一 writer
            micFrame = hw.reader();
            if isempty(micFrame)
                micFrame = zeros(blockSize, cfg.micNumChannels);
            end
            if size(micFrame,1) < blockSize && cfg.writeBlockPad
                micFrame(end+1:blockSize, :) = 0;
            end
            recorded(wrPtr:wrPtr+blockSize-1, :) = micFrame;
            wrPtr = wrPtr + blockSize;
        end
        recorded = recorded(1:totalSamples, :);

        if (cfg.saveFirstRaw && spk==1 && rep==1) || cfg.saveAllRaw
            rawRecordings{spk,rep} = recorded;
        end

        % 延迟估计（使用第一个麦克风）
        [cLag, lagsLag] = xcorr(recorded(:,1), sweepSig);
        [~, iLag] = max(cLag);
        delayCorr = lagsLag(iLag);
        repDelaysRaw(rep) = delayCorr;
        if rep == 1
            refDelay = delayCorr;
        end

        % 重复对齐（可选）
        if cfg.enableRepeatAlignment
            shiftSamples = delayCorr - refDelay;
            if abs(shiftSamples) > 2000
                fprintf('[warn]  Spk%d Rep%d 延迟差过大 shift=%d 标记剔除重复\n', spk, rep, shiftSamples);
                repDiscardRepeat(rep) = true;
            end
            if shiftSamples > 0
                recordedShift = [recorded(shiftSamples+1:end,:); zeros(shiftSamples, cfg.micNumChannels)];
            elseif shiftSamples < 0
                s = -shiftSamples;
                recordedShift = [zeros(s, cfg.micNumChannels); recorded(1:end-s,:)];
            else
                recordedShift = recorded;
            end
            [c2, l2] = xcorr(recordedShift(:,1), sweepSig);
            [~, i2] = max(c2);
            repDelaysCorr(rep) = l2(i2);
        else
            recordedShift = recorded;
            repDelaysCorr(rep) = delayCorr;
        end

        % 反卷积 & IR 提取
        irCurrent = zeros(Lh, cfg.micNumChannels);
        peakIdxEach = zeros(cfg.micNumChannels, 1);
        peakRelEach = false(cfg.micNumChannels, 1);
        snrMic = zeros(cfg.micNumChannels, 1);

        for m = 1:cfg.micNumChannels
            rec_m = recordedShift(:, m);
            outStruct = deconvolve_sweep(rec_m, sweepSig, cfg.fs, deconvParams);
            h_full = outStruct.h;
            pk = outStruct.peakIdx;
            reliable = isfield(outStruct, 'peakReliability') && outStruct.peakReliability;

            % 强制最小物理延迟
            if pk < cfg.minPhysDelaySamples
                searchStart = min(cfg.minPhysDelaySamples, length(h_full));
                searchEnd = min(searchStart + 400, length(h_full));
                if searchEnd > searchStart
                    [~, altRel] = max(abs(h_full(searchStart:searchEnd)));
                    pkAlt = searchStart + altRel - 1;
                    if pkAlt > pk && abs(h_full(pkAlt)) > 0.5 * max(abs(h_full))
                        pk = pkAlt;
                        reliable = true;
                        fprintf('[adjust] Spk%d Rep%d Mic%d 延迟重定位 %d -> %d\n', ...
                            spk, rep, m, outStruct.peakIdx, pk);
                    end
                end
            end
            if pk < cfg.minPhysDelaySamples
                reliable = false;
            end

            peakIdxEach(m) = pk;
            peakRelEach(m) = reliable;
            snrMic(m) = outStruct.snrEst;

            Luse = min(Lh, length(h_full));
            irCurrent(1:Luse, m) = h_full(1:Luse);
        end

        repPeakLists{rep} = peakIdxEach;
        repPeakReliable{rep} = peakRelEach;
        repSNRlist(rep, :) = snrMic;

        fprintf('[measure]  Spk%d Rep%d rawDelay=%d alignedDelay=%d 峰均值=%.1f reliableRatio=%.2f SNR(med)=%.2f discard=%d\n', ...
            spk, rep, delayCorr, repDelaysCorr(rep), mean(peakIdxEach), mean(peakRelEach), median(snrMic), repDiscardRepeat(rep));

        irRepList(:, :, rep) = irCurrent;
    end

    driftRaw = max(repDelaysRaw) - min(repDelaysRaw);
    driftCorr = max(repDelaysCorr) - min(repDelaysCorr);
    driftStable = driftCorr <= cfg.maxAllowedDriftSamples;
    if ~driftStable
        fprintf('[warn]  Spk%d 漂移不稳定 raw=%d corr=%d\n', spk, driftRaw, driftCorr);
    else
        fprintf('[measure] Spk%d 漂移稳定 raw=%d corr=%d\n', spk, driftRaw, driftCorr);
    end

    % 可靠峰值统计
    reliablePeaksAll = [];
    for rep = 1:cfg.repetitions
        if repDiscardRepeat(rep), continue; end
        pkList = repPeakLists{rep};
        relMask = repPeakReliable{rep};
        reliablePeaksAll = [reliablePeaksAll; pkList(relMask)];
    end

    reliableOutlierMask = false(size(reliablePeaksAll));
    reliableMedian = NaN; reliableIQR = NaN;
    if ~isempty(reliablePeaksAll)
        reliableMedian = median(reliablePeaksAll);
        q1 = prctile(reliablePeaksAll, 25);
        q3 = prctile(reliablePeaksAll, 75);
        reliableIQR = q3 - q1;
        if cfg.enableOutlierReject && reliableIQR > 0
            lowerFence = q1 - cfg.outlierIQRMultiplier * reliableIQR;
            upperFence = q3 + cfg.outlierIQRMultiplier * reliableIQR;
            reliableOutlierMask = reliablePeaksAll < lowerFence | reliablePeaksAll > upperFence;
            if any(reliableOutlierMask)
                fprintf('[filter] Spk%d 剔除 %d 个离群可靠峰\n', spk, sum(reliableOutlierMask));
            end
            reliablePeaksAll = reliablePeaksAll(~reliableOutlierMask);
            if ~isempty(reliablePeaksAll)
                reliableMedian = median(reliablePeaksAll);
                q1 = prctile(reliablePeaksAll, 25);
                q3 = prctile(reliablePeaksAll, 75);
                reliableIQR = q3 - q1;
            else
                reliableMedian = NaN;
                reliableIQR = NaN;
            end
        end
    end

    reliableCount = numel(reliablePeaksAll);
    if isempty(reliablePeaksAll)
        reliableMin = 0; reliableMax = 0;
    else
        reliableMin = min(reliablePeaksAll);
        reliableMax = max(reliablePeaksAll);
    end
    reliableRatioAll = sum(cellfun(@(x) sum(x), repPeakReliable)) / (cfg.repetitions * cfg.micNumChannels);

    % 平均 IR（原始）
    irAvgStruct = average_ir(irRepList, cfg.doAlignRepeats);
    irAvg = irAvgStruct.irAvg;

    % 每麦克风平均峰值
    peakPosAvgPerMic = zeros(cfg.micNumChannels, 1);
    for m = 1:cfg.micNumChannels
        [~, pkAvg] = max(abs(irAvg(:, m)));
        peakPosAvgPerMic(m) = pkAvg;
    end
    recommendedDelayFromAvg = round(median(peakPosAvgPerMic));

    if isempty(reliablePeaksAll)
        recommendedDelayReliable = recommendedDelayFromAvg;
        reliableMedian = recommendedDelayReliable;
        reliableIQR = NaN;
    else
        recommendedDelayReliable = round(reliableMedian);
    end

    if isempty(reliablePeaksAll)
        maxReliablePeak = recommendedDelayReliable;
    else
        maxReliablePeak = max(reliablePeaksAll);
    end
    recommendedFilterLen = min(cfg.irMaxLen, maxReliablePeak + floor(cfg.filterTailFraction * cfg.irMaxLen));

    % 尾部噪声评估
    tailLenEval = min(cfg.tailNoiseLen, size(irAvg,1));
    tailSeg = irAvg(end-tailLenEval+1:end, :);
    tailRMS = sqrt(mean(tailSeg(:).^2));

    % 导出对齐 IR（诊断用）
    if cfg.exportAlignedIR
        irAlignedStruct = average_ir(irRepList, true);
        irAligned = irAlignedStruct.irAvg;
        impulseResponsesAligned(:, :, spk) = irAligned;
    end

    % 可用性判定
    usable = driftStable && ...
             (reliableRatioAll >= cfg.reliableMinRatio) && ...
             (~isempty(reliablePeaksAll)) && ...
             ((isnan(reliableIQR) && reliableCount > 0) || reliableIQR <= cfg.reliableMaxIQR);

    fprintf('[summary-spk] Spk%d reliableCount=%d medianRel=%d IQR=%s min=%d max=%d ratio=%.2f avgDelay=%d relDelay=%d filtLen=%d usable=%d tailRMS=%.3e\n', ...
        spk, reliableCount, reliableMedian, num2str(reliableIQR), reliableMin, reliableMax, reliableRatioAll, ...
        recommendedDelayFromAvg, recommendedDelayReliable, recommendedFilterLen, usable, tailRMS);

    impulseResponses(:, :, spk) = irAvg;
    allIrRepList{spk} = irRepList;

    % 元数据
    spkMeta = struct();
    spkMeta.repDelaysRaw            = repDelaysRaw;
    spkMeta.repDelaysCorr           = repDelaysCorr;
    spkMeta.driftRaw                = driftRaw;
    spkMeta.driftCorr               = driftCorr;
    spkMeta.driftStable             = driftStable;
    spkMeta.irRepPeaks              = repPeakLists;
    spkMeta.irRepPeakReliable       = repPeakReliable;
    spkMeta.irRepSNR                = repSNRlist;
    spkMeta.irAvgPeakPos            = peakPosAvgPerMic;
    spkMeta.tailRMS                 = tailRMS;
    spkMeta.reliablePeaksAll        = reliablePeaksAll;
    spkMeta.reliableCount           = reliableCount;
    spkMeta.reliableMedian          = reliableMedian;
    spkMeta.reliableIQR             = reliableIQR;
    spkMeta.reliableMin             = reliableMin;
    spkMeta.reliableMax             = reliableMax;
    spkMeta.reliableRatioAll        = reliableRatioAll;
    spkMeta.recommendedDelayFromAvg = recommendedDelayFromAvg;
    spkMeta.recommendedDelayReliable= recommendedDelayReliable;
    spkMeta.recommendedFilterLen    = recommendedFilterLen;
    spkMeta.usable                  = usable;
    spkMeta.numReps                 = cfg.repetitions;

    meta.perSpeaker{spk} = spkMeta;
end

%% ============== 构建最终结构 ==============
secondary = struct();
secondary.impulseResponses = impulseResponses;
if cfg.exportAlignedIR
    secondary.impulseResponses_aligned = impulseResponsesAligned;
end
secondary.irRepList        = allIrRepList;
secondary.meta             = meta;
secondary.fs               = cfg.fs;
secondary.numSpeakers      = cfg.numSpeakers;
secondary.numMics          = cfg.micNumChannels;
secondary.irLength         = size(impulseResponses,1);
secondary.description      = 'Secondary path measurement v2.4 (ASIO4ALL single-device)';
secondary.measureConfig    = cfg;
secondary.sweepInfo        = sweepInfo;
secondary.deconvParams     = deconvParams;
secondary.timestampUtc     = datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS');

% 推荐延迟向量
recommendedVector = zeros(1, cfg.numSpeakers);
for spk = 1:cfg.numSpeakers
    recommendedVector(spk) = meta.perSpeaker{spk}.recommendedDelayReliable;
end
secondary.meta.recommendedDelaysVector = recommendedVector;
secondary.delayEstimateSamples = recommendedVector;

% 保存原始录音（可选）
if cfg.saveAllRaw || cfg.saveFirstRaw
    secondary.rawRecordingsInfo = struct('savedFirstRaw', cfg.saveFirstRaw, 'savedAllRaw', cfg.saveAllRaw);
    if cfg.saveAllRaw
        secondary.rawRecordings = rawRecordings;
    elseif cfg.saveFirstRaw
        secondary.rawRecordings = rawRecordings(1,1);
    end
end

% 保存文件
save(cfg.secondaryPathFile, 'secondary', '-v7.3');
fprintf('[measure] 保存完成 -> %s\n', cfg.secondaryPathFile);
fprintf('[measure] 推荐延迟向量: %s\n', mat2str(recommendedVector));

%% 汇总输出
for spk = 1:cfg.numSpeakers
    s = meta.perSpeaker{spk};
    fprintf('[final] Spk%d usable=%d relMedian=%d IQR=%s ratio=%.2f avgDelay=%d relDelay=%d filtLen=%d driftCorr=%d\n', ...
        spk, s.usable, s.reliableMedian, num2str(s.reliableIQR), s.reliableRatioAll, ...
        s.recommendedDelayFromAvg, s.recommendedDelayReliable, s.recommendedFilterLen, s.driftCorr);
end

try hw.release(); catch; end
fprintf('[measure] Done.\n');