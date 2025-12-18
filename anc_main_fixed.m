% 主仿真入口：多通道前馈 ANC + 在线次级路径建模 + 自适应延迟
clear; clc;
cfg = anc_config();
if ~check_anc_dependencies(cfg), return; end

%% ----------------- 加载输入音频 -----------------
[inputAudio, fs_in] = audioread(cfg.inputAudioFile);
if fs_in ~= cfg.fs
    error('输入采样率(%d)与配置(%d)不一致，请统一。', fs_in, cfg.fs);
end

% 确保通道数匹配
numRef = numel(cfg.micChannels.reference);
numErr = numel(cfg.micChannels.error);
expectedChans = numRef + numErr;
if size(inputAudio,2) ~= expectedChans
    error('输入音频通道数(%d) ≠ 预期(%d)', size(inputAudio,2), expectedChans);
end

totalSamples = size(inputAudio,1);
totalFrames = ceil(totalSamples / cfg.frameSize);
audioPadded = [inputAudio; zeros(totalFrames*cfg.frameSize - totalSamples, size(inputAudio,2))];

fprintf('✅ 输入加载完成: %d 样本, %d 帧, %d 通道\n', totalSamples, totalFrames, size(inputAudio,2));

%% ----------------- 加载反馈路径 F -----------------
loaded = load(cfg.feedbackPathFile, 'F');
if ~isfield(loaded, 'F')
    error('反馈路径文件 "%s" 未包含变量 F', cfg.feedbackPathFile);
end
F = loaded.F;
Lfb = size(F,1);
if ndims(F) ~= 3 || size(F,2) ~= numRef || size(F,3) ~= cfg.numSpeakers
    error('反馈路径 F 维度错误，应为 [Lfb, %d, %d]', numRef, cfg.numSpeakers);
end
fprintf('✅ 反馈路径 F 尺寸: [%d x %d x %d]\n', size(F));

%% ----------------- 加载次级路径 S -----------------
secData = load(cfg.secondaryPathFile);
if ~isfield(secData, 'secondary')
    error('次级路径文件缺少 "secondary" 结构体');
end
sec = secData.secondary;
if ~isfield(sec, 'impulseResponses_aligned')
    error('次级路径缺少 impulseResponses_aligned');
end
allIR = sec.impulseResponses_aligned;
errIdx = cfg.micChannels.error;
if size(allIR,2) < max(errIdx)
    error('次级路径麦克风数量不足');
end
S = allIR(:, errIdx, :);
[Ls, ne, ns] = size(S);
if ne ~= numErr || ns ~= cfg.numSpeakers
    error('次级路径维度不匹配: 期望 [%d, %d, %d], 实际 [%d, %d, %d]', Ls, numErr, cfg.numSpeakers, Ls, ne, ns);
end
fprintf('✅ 初始次级路径 S 尺寸: [%d x %d x %d]\n', size(S));

%% ----------------- 延迟估计初值 -----------------
if isfield(sec, 'delayEstimateSamples')
    delayD = sec.delayEstimateSamples + cfg.delayMarginSamples;
else
    delayD = round(0.0012 * cfg.fs) + cfg.delayMarginSamples;
end
delayD = max(delayD, cfg.minDelaySamples);
delayD = min(delayD, cfg.maxDelaySamples);
fprintf('ℹ️ 初始延迟 D=%d 样本 (≈ %.2f ms)\n', delayD, 1000*delayD/cfg.fs);

%% ----------------- 初始化 FXLMS 状态 -----------------
st = fxlms_recursive_init(cfg, S, delayD, numRef, numErr, cfg.numSpeakers);

% 反馈路径状态初始化（anc_main 里维护）
st.zi_fb = cell(numRef, cfg.numSpeakers);
for r = 1:numRef
    for s = 1:cfg.numSpeakers
        st.zi_fb{r,s} = zeros(max(Lfb-1,0), 1);
    end
end

% 初始化 xfHistory（确保存在并形状正确）
if cfg.timeFilterLen > 1
    Lw = cfg.timeFilterLen;
    st.xfHistory = zeros(Lw-1, numRef, cfg.numSpeakers, numErr);
end

% adaptive delay state placeholder
st.adaptiveDelayState = [];

%% ----------------- 滤波器系数预计算 -----------------
[b_bp,a_bp] = butter(4, [cfg.bandFreqLow cfg.bandFreqHigh]/(cfg.fs/2), 'bandpass');
bpStateLen = max(length(a_bp),length(b_bp))-1;
bpRefState = zeros(bpStateLen, numRef);
bpErrState = zeros(bpStateLen, numErr);

use_lowpass = isfield(cfg, 'ancLowpassHz') && ~isempty(cfg.ancLowpassHz);
if use_lowpass
    [b_lp, a_lp] = butter(4, cfg.ancLowpassHz / (cfg.fs/2));
end

%% ----------------- 日志预分配 -----------------
logFields = {'frame','refRms','errRms','improv','Wnorm','Wmax','mu','adaptEnable','outRms','clipping','delayD'};
logData = struct();
for f = 1:length(logFields)
    logData.(logFields{f}) = zeros(totalFrames,1);
end
logPtr = 1;

%% ----------------- 运行时变量 -----------------
frameCount = 0;
baselineCollected = false;
if ~isfield(cfg,'noiseBaselineFrames'), cfg.noiseBaselineFrames = 200; end
baselineBuffer = zeros(cfg.noiseBaselineFrames,1);
baselineCount = 0;
noiseBaselineRms = NaN;
outputBuffer = zeros(totalSamples, cfg.numSpeakers);

fprintf('【仿真启动】总帧=%d, 帧长=%d\n', totalFrames, cfg.frameSize);

%% ----------------- 主循环 -----------------
for frameIdx = 1:totalFrames
    startIdx = (frameIdx-1)*cfg.frameSize + 1;
    endIdx = startIdx + cfg.frameSize - 1;
    frameIn = audioPadded(startIdx:endIdx, :);
    
    refRaw = frameIn(:, cfg.micChannels.reference);
    errRaw = frameIn(:, cfg.micChannels.error);
    
    % 带通监控
    [refBand, bpRefState] = filter(b_bp,a_bp,refRaw,bpRefState);
    [errBand, bpErrState] = filter(b_bp,a_bp,errRaw,bpErrState);
    refRms = sqrt(mean(refBand(:).^2)+1e-12);
    errRms = sqrt(mean(errBand(:).^2)+1e-12);
    
    % 基线估计
    if ~baselineCollected
        baselineCount = baselineCount + 1;
        if baselineCount <= cfg.noiseBaselineFrames
            baselineBuffer(baselineCount) = errRms;
        end
        if frameCount >= cfg.noiseBaselineFrames && baselineCount >= cfg.noiseBaselineMinCount
            noiseBaselineRms = median(baselineBuffer(1:baselineCount));
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
    
    % 步长调度（简化）
    if frameCount <= cfg.warmupFrames
        current_mu = min(cfg.muInit + (cfg.muMax - cfg.muInit)*(frameCount/cfg.warmupFrames), cfg.muMax);
    else
        current_mu = cfg.muMax;
    end
    st.mu = current_mu;
    
    % 波束/参考使用
    if cfg.enableBeamforming
        mic_pos = (0:numRef-1) * cfg.refArraySpacing_m;
        ref_bf = reference_beamformer(refRaw, cfg.fs, mic_pos, 0);
        refUsed = repmat(ref_bf, 1, numRef);
    else
        refUsed = refRaw;
    end
    
    % 低通（若启用）
    if use_lowpass
        refUsed = filter(b_lp, a_lp, refUsed);
    end
    
    % 反馈泄漏抵消（若启用）
    if cfg.useFeedbackLeakCancel
        fbEst = zeros(size(refUsed));
        if frameIdx > 1
            lastOutStart = max(1, startIdx-cfg.frameSize);
            lastOut = outputBuffer(lastOutStart:startIdx-1, :);
            if size(lastOut,1) < cfg.frameSize
                lastOut = [zeros(cfg.frameSize-size(lastOut,1), cfg.numSpeakers); lastOut];
            end
        else
            lastOut = zeros(cfg.frameSize, cfg.numSpeakers);
        end
        for s = 1:cfg.numSpeakers
            for r = 1:numRef
                [fb_sig, st.zi_fb{r,s}] = filter(F(:,r,s), 1, lastOut(:,s), st.zi_fb{r,s});
                fbEst(:,r) = fbEst(:,r) + fb_sig;
            end
        end
        refUsed = refUsed - fbEst;
    end
    
    % 探测信号（按需生成）
    if cfg.enableOnlineSecPath
        probe_sig = cfg.probeAmplitude * randn(cfg.frameSize, cfg.numSpeakers);
    else
        probe_sig = zeros(cfg.frameSize, cfg.numSpeakers);
    end
    
    % FXLMS 处理（向量化函数）
    [y_out, st, dbg] = fxlms_recursive_process(st, refUsed, errRaw, cfg);
    y_total = y_out + probe_sig;
    
    % 输出控制 & 限幅
    if ~st.adaptEnable, y_total = y_total * 0.4; end
    y_total = y_total * outputRamp;
    peakPre = max(abs(y_total), [], 'all');
    clipping = false;
    if peakPre > cfg.maxOutput
        y_total = y_total * (cfg.maxOutput / peakPre);
        clipping = true;
    end
    y_total = max(min(y_total, cfg.maxOutput), -cfg.maxOutput);
    
    % 保存到输出缓冲
    actualEnd = min(endIdx, totalSamples);
    if startIdx <= totalSamples
        outputBuffer(startIdx:actualEnd, :) = y_total(1:(actualEnd-startIdx+1), :);
    end
    
    % 在线次级路径更新（按原逻辑/改进后的 estimator）
    if cfg.enableOnlineSecPath && mod(frameCount, cfg.secPathUpdateFrames)==0
        [S_new, st_sp] = online_sec_path_estimator([], probe_sig, errRaw, cfg);
        relChange = norm(S_new(:) - st.S(:)) / (norm(st.S(:)) + 1e-12);
        if relChange > getfieldOrDefault(cfg,'secPathUpdateRelChangeThreshold',0.5)
            % 安全更新：替换并重置相关状态以避免尺寸不匹配
            st.S = S_new;
            Ls_new = size(S_new,1);
            % 更新 st.zi 与 st.xfHistory 大小
            for r=1:numRef
                for s=1:cfg.numSpeakers
                    for e=1:numErr
                        st.zi{r,s,e} = zeros(max(0, Ls_new-1),1);
                        if cfg.timeFilterLen > 1
                            st.xfHistory(:,:,s,e) = zeros(cfg.timeFilterLen-1, numRef);
                        end
                    end
                end
            end
            fprintf('【在线S更新】帧=%d relChange=%.3f\n', frameCount, relChange);
        end
    end
    
    % ----------------- 自适应延迟估计集成 -----------------
    if cfg.enableAdaptiveDelay && frameCount > cfg.adaptiveDelayWarmupFrames && mod(frameCount, cfg.adaptiveDelayTriggerFrames) == 0
        [newD, st.adaptiveDelayState, ad_info] = adaptive_delay_estimator(st.adaptiveDelayState, refUsed, errRaw, cfg);
        % 如果 accepted 且与当前 delayD 有显著不同，则应用
        if isfield(ad_info,'accepted') && ad_info.accepted && newD ~= st.delayD
            prevD = st.delayD;
            delta = newD - prevD;
            % 限制单次最大变化（已在 estimator 平滑，但这里再次检查）
            maxChange = cfg.adaptiveDelayMaxChangeSamples;
            if abs(delta) > maxChange
                newD = prevD + sign(delta) * maxChange;
                delta = newD - prevD;
            end
            % 应用新延迟：为安全起见重置/对齐部分状态以避免跨帧不一致
            st.delayD = newD;
            % 重置参考缓冲与历史，以避免与新延迟不一致导致发散
            bufferLen = st.delayD + cfg.timeFilterLen;
            st.refBuffer = zeros(bufferLen, numRef);
            st.refIdx = 1;
            % 清空 xfHistory 与滤波器状态（保守策略）
            if cfg.timeFilterLen > 1
                st.xfHistory = zeros(cfg.timeFilterLen-1, numRef, cfg.numSpeakers, numErr);
            end
            for r = 1:numRef
                for s = 1:cfg.numSpeakers
                    for e = 1:numErr
                        st.zi{r,s,e} = zeros(max(0, length(st.S{e,s})-1),1);
                    end
                end
            end
            fprintf('【自适应延迟】帧=%d 接受新延迟: %d (prev=%d) conf=%.2f\n', frameCount, newD, prevD, ad_info.confidence);
            % 可选：短期内降低步长以缓冲转变带来的不稳定
            st.mu = st.mu * 0.5;
        else
            % 未采纳：可记录置信度
            if cfg.logging && isfield(ad_info,'confidence')
                % 将置信度写入日志（简单保存到 delayD 字段的负值做示例）
                % 更稳妥的做法是增加专门的日志结构记录 ad_info
            end
        end
    end
    
    % 改善量
    if baselineCollected
        improv_dB = 20*log10(noiseBaselineRms/(errRms+1e-12));
    else
        improv_dB = NaN;
    end
    
    % 日志记录
    if cfg.logging
        logData.frame(logPtr) = frameCount;
        logData.refRms(logPtr) = refRms;
        logData.errRms(logPtr) = errRms;
        logData.improv(logPtr) = improv_dB;
        logData.Wnorm(logPtr) = dbg.W_norm;
        logData.Wmax(logPtr) = dbg.W_maxAbs;
        logData.mu(logPtr) = st.mu;
        logData.adaptEnable(logPtr) = st.adaptEnable;
        logData.outRms(logPtr) = sqrt(mean(y_total(:).^2)+1e-12);
        logData.clipping(logPtr) = clipping;
        logData.delayD(logPtr) = st.delayD;
        logPtr = logPtr + 1;
    end
    
    % 状态打印
    if mod(frameCount, cfg.statusPrintFrames) == 0
        phaseStr = {'静音','斜坡','正常'};
        ps = phaseStr{min(adaptPhase+1,3)};
        if clipping, clipStr = '[限幅]'; else clipStr = ''; end
        fprintf('【帧=%4d】阶段=%s | Ref=%.2e Err=%.2e 改善=%.1f dB | μ=%.1e | ‖W‖=%.3e | delayD=%d %s\n', ...
            frameCount, ps, refRms, errRms, improv_dB, st.mu, dbg.W_norm, st.delayD, clipStr);
    end
end

%% ----------------- 结果汇总与保存 -----------------
if cfg.logging
    recordedN = max(0, logPtr-1);
else
    recordedN = 0;
end

if baselineCollected && recordedN>0
    lastN = min(100, recordedN);
    finalErrAvg = mean(logData.errRms(recordedN-lastN+1:recordedN));
    finalImprov = 20*log10(noiseBaselineRms/(finalErrAvg+1e-12));
    validImprov = logData.improv(1:recordedN);
    validImprov = validImprov(~isnan(validImprov) & validImprov>-50);
    fprintf('\n【最终结果】基线=%.3e | 最终误差=%.3e | 改善=%.2f dB | 平均改善=%.2f dB | 最大改善=%.2f dB\n', ...
        noiseBaselineRms, finalErrAvg, finalImprov, mean(validImprov), max(validImprov));
else
    fprintf('\n⚠️ 未成功锁定噪声基线或无日志数据，无法计算最终改善。\n');
end

if cfg.logging
    trimmedLog = struct();
    for f = 1:length(logFields)
        fld = logFields{f};
        trimmedLog.(fld) = logData.(fld)(1:recordedN);
    end
    saveStruct = struct('logData', trimmedLog, 'cfg', cfg, 'outputBuffer', outputBuffer);
    save(cfg.logFile, 'saveStruct');
    fprintf('✅ 日志保存: %s\n', cfg.logFile);
    
    outAudio = outputBuffer;
    sc = max(abs(outAudio(:)));
    if sc>0, outAudio = outAudio/sc*0.95; end
    audiowrite(cfg.outputAudioFile, outAudio, cfg.fs);
    fprintf('✅ 输出音频保存: %s\n', cfg.outputAudioFile);
end

if cfg.logging
    anc_plot_results(trimmedLog);
end

fprintf('\n【仿真完成】\n');

%% ----------------- 辅助函数 -----------------
function cfg = validate_fields(cfg, fields, defaultValue)
    for i = 1:length(fields)
        if ~isfield(cfg, fields{i})
            cfg.(fields{i}) = defaultValue;
        end
    end
end

function val = getfieldOrDefault(s, field, default)
    if isstruct(s) && isfield(s, field), val = s.(field); else val = default; end
end