% record_secondary_path.m 多扬声器→多误差麦克风次级路径测量（ESS）
% 版本：v3.0 彻底修复版 - 解决互相关失败和延迟估计问题

clear; clc;

%% ============== 配置区 ==============
cfg = anc_config();

% 强制检查并修正配置
if cfg.minPhysDelaySamples < 50
    cfg.minPhysDelaySamples = 50;
    fprintf('[config] 警告: minPhysDelaySamples调整到50\n');
end
if cfg.maxPhysDelaySamples < cfg.minPhysDelaySamples + 1000
    cfg.maxPhysDelaySamples = cfg.minPhysDelaySamples + 10000;
    fprintf('[config] 警告: maxPhysDelaySamples调整到%d\n', cfg.maxPhysDelaySamples);
end

% 仅处理误差麦克风
errMicIdx = cfg.micChannels.error;
numErrMics = length(errMicIdx);
fprintf('[measure] 误差麦克风通道: %s\n', mat2str(errMicIdx));

%% ============== 初始化硬件 ==============
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
    sweepCoreLen = length(sweepDrive);
    fprintf('[measure] Sweep 总长=%d 样本，核心sweep=%d 样本\n', totalSamples, sweepCoreLen);

    %% ============== 容器初始化 ==============
    Lh = cfg.irMaxLen;
    impulseResponses = zeros(Lh, numErrMics, cfg.numSpeakers);
    impulseResponsesAligned = zeros(Lh, numErrMics, cfg.numSpeakers);

    blockSize = cfg.timeFrameSamples;
    numBlocks = ceil(totalSamples / blockSize);

    % 元数据初始化
    meta = struct();
    meta.global = struct('sampleRate',cfg.fs,'sweepLen',totalSamples,'repetitions',cfg.repetitions);
    meta.perSpeaker = cell(cfg.numSpeakers,1);

    %% ============== 测量循环 ==============
    hw = [];
    for spk = 1:cfg.numSpeakers
        fprintf('\n[measure] ===== 扬声器 %d =====\n', spk);

        repRecordedRaw = cell(cfg.repetitions,1);
        repPeakLists_raw = cell(cfg.repetitions,1);
        repPeakReliable = cell(cfg.repetitions,1);
        repSNRlist = zeros(cfg.repetitions, numErrMics);
        irRepList = zeros(Lh, numErrMics, cfg.repetitions);

        % === 采集所有 repeats ===
        for rep = 1:cfg.repetitions
            fprintf('[measure]  播放重复 %d/%d...\n', rep, cfg.repetitions);
            
            if rep == 1
                fprintf('[measure] 初始化硬件...\n');
                hw = hardware_init_measure(cfg);
            else
                hw.release();
                pause(0.3);
                hw = hardware_init_measure(cfg);
            end

            % 播放信号
            spkDriveSig = cfg.spkAmplitude(spk) * sweepSig(:);

            % 防削波
            maxVal = max(abs(spkDriveSig));
            if maxVal > 0.98
                scale = 0.95 / maxVal;
                spkDriveSig = spkDriveSig * scale;
                fprintf('[measure] 播放缩放 %.2f\n', scale);
            end

            % 初始化录音buffer
            ptr = 1;
            recordedFull = zeros(numBlocks * blockSize, cfg.micNumChannels);

            % 预热
            for pr = 1:cfg.preRollFrames
                hw.writer(zeros(blockSize, cfg.numSpeakers));
                hw.reader();
            end

            % 主播放循环
            rdPtr = 1;
            for b = 1:numBlocks
                outBlock = zeros(blockSize, cfg.numSpeakers);
                if ptr <= totalSamples
                    nToCopy = min(blockSize, totalSamples - ptr + 1);
                    outBlock(1:nToCopy, spk) = spkDriveSig(ptr : ptr + nToCopy - 1);
                    ptr = ptr + nToCopy;
                end
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

            % 确保录音长度正确
            if size(recordedFull, 1) ~= totalSamples
                fprintf('[WARN] Spk%d Rep%d: 录音长度不匹配 (期望=%d, 实际=%d)\n', ...
                    spk, rep, totalSamples, size(recordedFull, 1));
                
                if size(recordedFull, 1) > totalSamples
                    recordedFull = recordedFull(1:totalSamples, :);
                else
                    recordedFull(end+1:totalSamples, :) = 0;
                end
            end
            
            recorded = recordedFull(:, errMicIdx);
            repRecordedRaw{rep} = recorded;
        end
        hw.release();

        % === 关键修复：基于能量检测确定信号起始 ===
        % 不再依赖互相关，直接寻找能量峰值
        
        % 计算所有录音的平均能量
        allRec = cat(2, repRecordedRaw{:});
        avgRec = mean(allRec, 2);  % 平均所有repeats
        
        % 使用滑动窗口检测能量峰值
        winSize = 1000;
        energy = movmean(sum(avgRec.^2, 2), winSize);
        
        % 找到能量超过阈值10%的第一个点
        energy_threshold = 0.1 * max(energy);
        signal_start = find(energy > energy_threshold, 1, 'first');
        
        if isempty(signal_start)
            signal_start = cfg.minPhysDelaySamples;
            fprintf('[WARN] Spk%d: 未检测到明显信号，使用默认延迟\n', spk);
        else
            fprintf('[info] Spk%d: 能量检测到信号起始于样本 %d\n', spk, signal_start);
        end
        
        % 确定截取窗口
        extractStart = max(1, signal_start - 500);  % 提前500个样本
        extractEnd = min(size(avgRec,1), extractStart + sweepCoreLen + Lh*2 - 1);
        
        fprintf('[info] Spk%d: 截取窗口 [%d, %d]\n', spk, extractStart, extractEnd);

        for rep = 1:cfg.repetitions
            rec = repRecordedRaw{rep};
            
            % 截取窗口
            if extractEnd > size(rec,1)
                recExtract = [rec(extractStart:end, :); zeros(extractEnd-size(rec,1), size(rec,2))];
            else
                recExtract = rec(extractStart:extractEnd, :);
            end
            
            offsetInRec = extractStart - 1;

            % 对每个麦克风进行反卷积
            peakIdxEach_raw = zeros(numErrMics, 1);
            peakRelEach = false(numErrMics, 1);
            snrMic = zeros(numErrMics, 1);
            irCurrent = zeros(Lh, numErrMics);

            for m = 1:numErrMics
                % 使用改进的反卷积函数
                outStruct = deconvolve_sweep(recExtract(:,m), sweepDrive, cfg);
                h_full = outStruct.h;
                pk_local = outStruct.peakIdx;

                % 映射到全局索引
                pk_global = offsetInRec + pk_local;
                
                % 质量检查
                delay_ok = (pk_global >= cfg.minPhysDelaySamples) && ...
                           (pk_global <= cfg.maxPhysDelaySamples * 10);  % 放宽上限
                snr_ok = (outStruct.snrEst >= cfg.snrThresholdDB);
                reliable = delay_ok && snr_ok && outStruct.peakReliability;

                fprintf('  [INFO] Mic%d: pk_global=%d, SNR=%.2f dB, reliable=%d\n', ...
                    m, pk_global, outStruct.snrEst, reliable);

                peakIdxEach_raw(m) = pk_global;
                peakRelEach(m) = reliable;
                snrMic(m) = outStruct.snrEst;

                % 保存IR
                Luse = min(Lh, length(h_full));
                irCurrent(1:Luse, m) = h_full(1:Luse);
            end

            repPeakLists_raw{rep} = peakIdxEach_raw;
            repPeakReliable{rep} = peakRelEach;
            repSNRlist(rep, :) = snrMic;
            irRepList(:, :, rep) = irCurrent;
        end

        % === 平均IR ===
        irAvg = mean(irRepList, 3);
        impulseResponses(:, :, spk) = irAvg;

        % === 可用性判定 ===
        allPeaks = [];
        for rep = 1:cfg.repetitions
            peaks = repPeakLists_raw{rep}(repPeakReliable{rep});
            allPeaks = [allPeaks; peaks(:)];
        end

        if isempty(allPeaks)
            usable = false;
            recommendedDelay = cfg.minPhysDelaySamples + 100;
            fprintf('[diag] Spk%d: 没有可靠的峰值\n', spk);
        else
            % 计算统计指标
            peak_std = std(allPeaks);
            peak_range = range(allPeaks);
            median_peak = median(allPeaks);
            median_snr = median(repSNRlist(repSNRlist > -Inf));
            
            % 判定条件
            stable_peaks = (peak_std < 100) && (peak_range < 500);
            good_snr = (median_snr >= cfg.snrThresholdDB);
            usable = stable_peaks && good_snr;
            
            recommendedDelay = round(median_peak);
            
            fprintf('[diag] Spk%d: 峰值中位数=%d, 标准差=%.1f, 范围=%d, SNR中位数=%.2f dB, usable=%d\n', ...
                spk, median_peak, peak_std, peak_range, median_snr, usable);
        end

        % 保存元数据
        meta.perSpeaker{spk} = struct(...
            'repPeakPositions', repPeakLists_raw, ...
            'repPeakReliable', repPeakReliable, ...
            'repSNRlist', repSNRlist, ...
            'recommendedDelay', recommendedDelay, ...
            'usable', usable, ...
            'extractWindow', [extractStart, extractEnd] ...
            );
    end

    %% ============== 构建输出结构 ==============
    secondary = struct();
    secondary.impulseResponses = impulseResponses;
    secondary.fs = cfg.fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;
    secondary.errorMicPhysicalChannels = errMicIdx;
    secondary.irLength = Lh;
    secondary.description = 'Secondary path v3.0 (fixed correlation failure)';
    secondary.timestampUtc = datestr(datetime('now','TimeZone','UTC'),'yyyy-mm-dd HH:MM:SS');

    % 推荐延迟
    recommendedVector = zeros(1, cfg.numSpeakers);
    for i = 1:cfg.numSpeakers
        if meta.perSpeaker{i}.usable
            recommendedVector(i) = meta.perSpeaker{i}.recommendedDelay;
        else
            recommendedVector(i) = cfg.minPhysDelaySamples + 100;
        end
    end
    secondary.delayEstimateSamples = recommendedVector;
    
    % 保存诊断信息
    if cfg.saveDiagnosticInfo
        secondary.meta = meta;
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
