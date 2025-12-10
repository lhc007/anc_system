% 主仿真入口：多通道前馈 ANC + 在线次级路径建模（数组化 W 版）
clear; clc;

cfg = anc_config();
if ~check_anc_dependencies(cfg), return; end

%% 加载输入噪声（多通道）
[inputAudio, fs_in] = audioread(cfg.inputAudioFile);
if fs_in ~= cfg.fs
    error('输入采样率(%d)与配置(%d)不一致，请统一。', fs_in, cfg.fs);
end

expectedChans = numel(cfg.micChannels.reference) + numel(cfg.micChannels.error);
if size(inputAudio,2) ~= expectedChans
    error('输入音频通道数(%d) ≠ 预期(%d)', size(inputAudio,2), expectedChans);
end

totalSamples = size(inputAudio,1);
frameSize    = cfg.frameSize;
totalFrames  = ceil(totalSamples / frameSize);

audioPadded  = [inputAudio; zeros(totalFrames*frameSize - totalSamples, size(inputAudio,2))];
numRef = numel(cfg.micChannels.reference);
numErr = numel(cfg.micChannels.error);

fprintf('✅ 输入加载完成: %d 样本, %d 帧, %d 通道\n', totalSamples, totalFrames, size(inputAudio,2));

%% 加载反馈路径 F: [Lfb x numRef x numSpeakers]
loaded = load(cfg.feedbackPathFile, 'F');
if ~isfield(loaded, 'F'), error('反馈路径文件 "%s" 未包含变量 F', cfg.feedbackPathFile); end
F = loaded.F;
if ndims(F) ~= 3 || size(F,2) ~= numRef || size(F,3) ~= cfg.numSpeakers
    error('反馈路径 F 维度错误，应为 [Lfb, %d, %d]', numRef, cfg.numSpeakers);
end
Lfb = size(F,1);
fprintf('✅ 反馈路径 F 尺寸: [%d x %d x %d]\n', size(F));

%% 加载次级路径 S: [Ls x numErr x numSpeakers]
secData = load(cfg.secondaryPathFile);
if ~isfield(secData, 'secondary'), error('次级路径文件缺少 "secondary" 结构体'); end
sec = secData.secondary;
if ~isfield(sec, 'impulseResponses_aligned')
    error('次级路径缺少 impulseResponses_aligned');
end
allIR = sec.impulseResponses_aligned;
errIdx = cfg.micChannels.error;
if size(allIR,2) < max(errIdx), error('次级路径麦克风数量不足'); end
S = allIR(:, errIdx, :);
[Ls, ne, ns] = size(S);
if ne ~= numErr || ns ~= cfg.numSpeakers
    error('次级路径维度不匹配: 期望 [%d, %d, %d], 实际 [%d, %d, %d]', Ls, numErr, cfg.numSpeakers, Ls, ne, ns);
end
fprintf('✅ 初始次级路径 S 尺寸: [%d x %d x %d]\n', size(S));

%% 延迟估计
if isfield(sec, 'delayEstimateSamples')
    delayD = sec.delayEstimateSamples + cfg.delayMarginSamples;
else
    delayD = round(0.0012 * cfg.fs) + cfg.delayMarginSamples;
end
delayD = min(max(delayD, cfg.minDelaySamples), cfg.maxDelaySamples);
fprintf('ℹ️ 初始延迟 D=%d 样本 (≈ %.2f ms)\n', delayD, 1000*delayD/cfg.fs);

%% 初始化 FXLMS（数组化 W 版本）
st = fxlms_recursive_init(cfg, S, delayD, numRef, ne, ns);

% 初始化反馈路径状态
st.zi_fb = cell(numRef, cfg.numSpeakers);
for r = 1:numRef
    for s = 1:cfg.numSpeakers
        st.zi_fb{r,s} = zeros(Lfb-1, 1);
    end
end

%% 监控带通滤波器
[b_bp,a_bp] = butter(4, [cfg.bandFreqLow cfg.bandFreqHigh]/(cfg.fs/2),'bandpass');
bpRefState = zeros(max(length(a_bp),length(b_bp))-1, numRef);
bpErrState = zeros(max(length(a_bp),length(b_bp))-1, numErr);

%% 功能开关默认值保护
cfg = validate_fields(cfg, {'enableOnlineSecPath', 'enableAdaptiveDelay', 'enableBeamforming'}, false);

%% 预分配日志
logFields = {'frame','refRms','errRms','improv','Wnorm','Wmax','mu','adaptEnable','outRms','clipping','delayD'};
for f = 1:length(logFields)
    logData.(logFields{f}) = zeros(totalFrames, 1);
end
logPtr = 1;

%% 初始化运行时变量
frameCount = 0;
baselineCollected = false;
baselineBuffer = [];
noiseBaselineRms = NaN;
outputBuffer = zeros(totalSamples, cfg.numSpeakers);

fprintf('【仿真启动】总帧=%d, 帧长=%d\n', totalFrames, frameSize);

%% 主循环
for frameIdx = 1:totalFrames
    startIdx = (frameIdx-1)*frameSize + 1;
    endIdx   = startIdx + frameSize - 1;
    frameIn  = audioPadded(startIdx:endIdx, :);

    refRaw = frameIn(:, cfg.micChannels.reference);
    errRaw = frameIn(:, cfg.micChannels.error);

    % 带通 RMS 监控
    [refBand, bpRefState] = filter(b_bp,a_bp,refRaw,bpRefState);
    [errBand, bpErrState] = filter(b_bp,a_bp,errRaw,bpErrState);
    refRms = sqrt(mean(refBand(:).^2)+1e-12);
    errRms = sqrt(mean(errBand(:).^2)+1e-12);

    % 噪声基线估计
    if ~baselineCollected
        baselineBuffer(end+1) = errRms;
        if frameCount >= cfg.noeseBaselineFrames && numel(baselineBuffer) >= cfg.noiseBaselineMinCount
            noiseBaselineRms = median(baselineBuffer);
            baselineCollected = true;
            fprintf('✅ 噪声基线锁定: %.3e\n', noiseBaselineRms);
        end
    end

    frameCount = frameCount + 1;

    % 时序控制
    if frameCount <= cfg.startMuteFrames
        adaptPhase = 0; outputRamp = 0;
    elseif frameCount <= cfg.startMuteFrames + cfg.rampFrames
        adaptPhase = 1; outputRamp = (frameCount - cfg.startMuteFrames)/cfg.rampFrames;
    else
        adaptPhase = 2; outputRamp = 1;
    end
    st.adaptEnable = (adaptPhase >= 2);

    % 步长调度
    if frameCount <= cfg.warmupFrames
        current_mu = min(cfg.muInit + (cfg.muMax - cfg.muInit)*(frameCount/cfg.warmupFrames), cfg.muMax);
    else
        current_mu = cfg.muMax;
    end
    st.mu = current_mu;

    % 波束成形（可选）
    if cfg.enableBeamforming
        mic_pos = (0:numRef-1) * cfg.refArraySpacing_m;
        ref_bf = reference_beamformer(refRaw, cfg.fs, mic_pos, 0);
        refUsed = repmat(ref_bf, 1, numRef);
    else
        refUsed = refRaw;
    end

    % 反馈泄漏抵消
    if cfg.useFeedbackLeakCancel
        fbEst = zeros(size(refUsed));
        if frameIdx > 1
            lastOut = outputBuffer(max(1,startIdx-frameSize):startIdx-1, :);
            if size(lastOut,1) < frameSize
                lastOut = [zeros(frameSize-size(lastOut,1), cfg.numSpeakers); lastOut];
            end
        else
            lastOut = zeros(frameSize, cfg.numSpeakers);
        end
        for s = 1:cfg.numSpeakers
            for r = 1:numRef
                [fb_sig, st.zi_fb{r,s}] = filter(F(:,r,s), 1, lastOut(:,s), st.zi_fb{r,s});
                fbEst(:,r) = fbEst(:,r) + fb_sig;
            end
        end
        refUsed = refUsed - fbEst;
    end

    % 探测信号
    probe_sig = cfg.enableOnlineSecPath * cfg.probeAmplitude * randn(frameSize, cfg.numSpeakers);

    % FXLMS 处理（数组化 W + 窗口向量化）
    [y_out, st, dbg] = fxlms_recursive_process(st, refUsed, errRaw, cfg);
    y_total = y_out + probe_sig;

    % 输出控制
    if ~st.adaptEnable, y_total = y_total * 0.4; end
    y_total = y_total * outputRamp;
    peakPre = max(abs(y_total), [], 'all');
    clipping = false;
    if peakPre > cfg.maxOutput
        y_total = y_total * (cfg.maxOutput / peakPre);
        clipping = true;
    end
    y_total = max(min(y_total, cfg.maxOutput), -cfg.maxOutput);

    % 保存有效输出
    actualEnd = min(endIdx, totalSamples);
    if startIdx <= totalSamples
        outputBuffer(startIdx:actualEnd, :) = y_total(1:(actualEnd-startIdx+1), :);
    end

    % 在线次级路径更新（示意）
    if cfg.enableOnlineSecPath && mod(frameCount, 100) == 0
        [S_new, ~] = online_sec_path_estimator([], probe_sig, errRaw, cfg);
        % 简单一致性检查（相对变化不大才更新）
        if norm(S_new(:) - S(:)) / (norm(S(:)) + 1e-12) < 0.5
            S = S_new; % 更新原始 S
            % 同步 st.S_cell 与 zi_fxf
            for e = 1:numErr
                for s = 1:cfg.numSpeakers
                    ir_raw = squeeze(S(:, e, s));
                    lastNonZero = find(ir_raw ~= 0, 1, 'last');
                    if ~isempty(lastNonZero), ir = ir_raw(1:lastNonZero); else, ir = 0; end
                    st.S_cell{e, s} = ir(:);
                    Ls_new = length(ir);
                    for r = 1:numRef
                        st.zi_fxf{e, s, r} = zeros(max(0, Ls_new - 1), 1);
                    end
                end
            end
            fprintf('【在线S更新】帧=%d\n', frameCount);
        end
    end

    % 自适应延迟（略）

    % 改善量
    if baselineCollected
        improv_dB = 20*log10(noiseBaselineRms/(errRms+1e-12));
    else
        improv_dB = NaN;
    end

    % 日志记录
    if cfg.logging
        logData.frame(logPtr)       = frameCount;
        logData.refRms(logPtr)      = refRms;
        logData.errRms(logPtr)      = errRms;
        logData.improv(logPtr)      = improv_dB;
        logData.Wnorm(logPtr)       = dbg.W_norm;
        logData.Wmax(logPtr)        = dbg.W_maxAbs;
        logData.mu(logPtr)          = st.mu;
        logData.adaptEnable(logPtr) = st.adaptEnable;
        logData.outRms(logPtr)      = sqrt(mean(y_total(:).^2)+1e-12);
        logData.clipping(logPtr)    = clipping;
        logData.delayD(logPtr)      = st.delayD;
        logPtr = logPtr + 1;
    end

    % 状态打印
    if mod(frameCount, cfg.statusPrintFrames)==0
        phaseStr = {'静音','斜坡','正常'};
        ps = phaseStr{min(adaptPhase+1,3)};
        clipStr = ''; if clipping, clipStr = '[限幅]'; end
        outRms = sqrt(mean(y_total(:).^2)+1e-12);
        fprintf('【帧=%4d】阶段=%s | Ref=%.2e Err=%.2e 改善=%.1f dB | μ=%.1e | ‖W‖=%.3e | delayD=%d %s\n', ...
            frameCount, ps, refRms, errRms, improv_dB, st.mu, dbg.W_norm, st.delayD, clipStr);
    end
end

%% 结果汇总
if baselineCollected
    lastN = min(100, logPtr-1);
    finalErrAvg = mean(logData.errRms(logPtr-lastN:logPtr-1));
    finalImprov = 20*log10(noiseBaselineRms/(finalErrAvg+1e-12));
    validImprov = logData.improv(1:logPtr-1);
    validImprov = validImprov(~isnan(validImprov) & validImprov>-50);
    fprintf('\n【最终结果】基线=%.3e | 最终误差=%.3e | 改善=%.2f dB | 平均改善=%.2f dB | 最大改善=%.2f dB\n', ...
        noiseBaselineRms, finalErrAvg, finalImprov, mean(validImprov), max(validImprov));
else
    fprintf('\n⚠️ 未成功锁定噪声基线，无法计算最终改善。\n');
end

%% 保存结果
if cfg.logging
    saveStruct = struct('logData', logData(1:logPtr-1), 'cfg', cfg, 'outputBuffer', outputBuffer);
    save(cfg.logFile, 'saveStruct');
    fprintf('✅ 日志保存: %s\n', cfg.logFile);

    outAudio = outputBuffer;
    sc = max(abs(outAudio(:)));
    if sc>0, outAudio = outAudio/sc*0.95; end
    audiowrite(cfg.outputAudioFile, outAudio, cfg.fs);
    fprintf('✅ 输出音频保存: %s\n', cfg.outputAudioFile);
end

anc_plot_results(logData(1:logPtr-1));
fprintf('\n【仿真完成】\n');

%% 辅助函数
function cfg = validate_fields(cfg, fields, defaultValue)
    for i = 1:length(fields)
        if ~isfield(cfg, fields{i})
            cfg.(fields{i}) = defaultValue;
        end
    end
end