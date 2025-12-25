function record_secondary_path_v3()
% record_secondary_path_v3.m 多通道管道ANC次级路径测量系统
% 版本：v4.2 优化版
% 特性：自适应延迟估计、多通道同步处理、实时质量监控、抗混叠处理
% 修正重点：
%   - 修复延迟估计中的边界条件错误
%   - 添加相干性检查功能
%   - 改进错误处理逻辑
%   - 移除IR幅度归一化（保留物理增益）

clear; clc;

%% ============== 配置加载与验证 ==============
cfg = anc_config();

% 配置验证与自动修正
fprintf('[init] 加载配置...\n');
cfg = validate_config(cfg);

% 误差麦克风通道
errMicIdx = cfg.micChannels.error;
numErrMics = length(errMicIdx);
fprintf('[init] 误差麦克风通道: %s\n', mat2str(errMicIdx));

%% ============== 初始化系统 ==============
try
    %% ============== 生成激励信号 ==============
    fprintf('[signal] 生成ESS激励信号...\n');
    [sweepCore, ~, ~] = generate_sweep_v3(cfg);
    fs = cfg.fs();
    
    % --- 第一层：为 ESS 反卷积添加必要的前后静音（来自 cfg.padLeading/Trailing）---
    padLeadSamples = round(cfg.padLeading * fs);
    padTrailSamples = round(cfg.padTrailing * fs);
    sweepWithEssPad = [zeros(padLeadSamples, 1); ...
                       sweepCore; ...
                       zeros(padTrailSamples, 1)];
    
    % --- 可选：低频增强（作用于纯扫频或带 pad 信号？建议作用于 sweepCore）---
    if cfg.enableLowFreqBoost
        % 示例：在 sweepCore 上叠加低频正弦（需确保长度一致）
        t_core = (0:length(sweepCore)-1)' / fs;
        lowFreqTone = 0.1 * sin(2*pi*20*t_core);  % 20Hz, 幅度0.1
        sweepCoreBoosted = sweepCore + lowFreqTone;
        
        % 重新构建带 pad 的信号
        sweepWithEssPad = [zeros(padLeadSamples, 1); ...
                           sweepCoreBoosted; ...
                           zeros(padTrailSamples, 1)];
    else
        % 已构建 sweepWithEssPad
    end
    
    % --- 归一化驱动信号 ---
    maxAmp = max(abs(sweepWithEssPad));
    if maxAmp > 0
        sweepDrive = sweepWithEssPad / maxAmp * cfg.amplitude;
    else
        sweepDrive = sweepWithEssPad * cfg.amplitude;
    end
    
    % --- 第二层：为测量添加额外的前后静音（preSilenceSec / postSilenceSec）---
    preSilence = zeros(round(cfg.preSilenceSec * fs), 1);
    postSilence = zeros(round(cfg.postSilenceSec * fs), 1);
    sweepSig = [preSilence; sweepDrive; postSilence];
    
    % ========== 日志打印 ==========
    sweepCoreLen = length(sweepCore);          % 真正的核心扫频长度
    essTotalLen = length(sweepDrive);          % ESS激励总长（含ESS-pad）
    totalSamples = length(sweepSig);           % 最终播放总长（含 pre/post silence）
    
    fprintf('[signal] 信号参数:\n');
    fprintf('  - 核心扫频时长: %.3f s (%d 样本)\n', cfg.sweepDuration, sweepCoreLen);
    fprintf('  - ESS激励长度（含内部pad）: %d 样本\n', essTotalLen);
    fprintf('  - 测量信号总长度: %d 样本\n', totalSamples);
    fprintf('  - 采样率: %d Hz\n', fs);

    %% ============== 容器初始化 ==============
    Lh = cfg.irMaxLen;
    impulseResponses = zeros(Lh, numErrMics, cfg.numSpeakers);
    
    % 数据结构
    meta = struct();
    meta.global = struct(...
        'sampleRate', cfg.fs, ...
        'sweepLen', totalSamples, ...
        'sweepCoreLen', sweepCoreLen, ...
        'repetitions', cfg.repetitions, ...
        'irLength', Lh, ...
        'errorMicChannels', errMicIdx, ...
        'timestamp', datetime('now', 'TimeZone', 'local'));
    
    meta.perSpeaker = cell(cfg.numSpeakers, 1);

    %% ============== 主测量循环 ==============
    fprintf('\n[measure] 开始次级路径测量...\n');
    
    for spk = 1:cfg.numSpeakers
        fprintf('\n%s\n', repmat('=', 1, 60));
        fprintf('[measure] 扬声器 %d/%d\n', spk, cfg.numSpeakers);
        fprintf('%s\n', repmat('=', 1, 60));
        
        % 扬声器特定数据结构
        spkData = struct();
        spkData.irRaw = cell(cfg.repetitions, 1);
        spkData.delayEst = zeros(cfg.repetitions, 1);
        spkData.snrEst = zeros(cfg.repetitions, numErrMics);
        spkData.peakPos = zeros(cfg.repetitions, numErrMics);
        spkData.coherenceEst = zeros(cfg.repetitions, numErrMics); % 新增相干性矩阵
        spkData.irFinal = zeros(Lh, numErrMics, cfg.repetitions);
        
        % === 阶段1: 数据采集 ===
        fprintf('[measure] 阶段1: 数据采集\n');
        
        for rep = 1:cfg.repetitions
            fprintf('  重复 %d/%d... ', rep, cfg.repetitions);
            
            % 硬件初始化
            if rep == 1
                hw = hardware_init_measure(cfg);
                fprintf('硬件已初始化\n');
            else
                pause(0.5); % 重连间隔
                hw = hardware_init_measure(cfg);
            end
            
            % 生成播放信号（带防削波）
            spkDriveSig = cfg.spkAmplitude(spk) * sweepSig(:);
            maxVal = max(abs(spkDriveSig));
            if maxVal > 0.95
                scale = 0.9 / maxVal;
                spkDriveSig = spkDriveSig * scale;
                fprintf('    防削波缩放: %.2f\n', scale);
            end
            
            % 播放并录音
            [recordedFull, playbackInfo] = play_and_record(hw, spkDriveSig, spk, cfg);
            
            % 提取误差麦克风通道
            recorded = recordedFull(:, errMicIdx);
            
            % 保存原始数据
            spkData.irRaw{rep} = recorded;
            
            % 释放硬件
            hw.release();
            
            fprintf('完成\n');
        end
        
        % === 阶段2: 自适应延迟估计 ===
        fprintf('[measure] 阶段2: 自适应延迟估计\n');
        
        % 使用所有重复的平均信号进行延迟估计
        allRecorded = zeros(size(spkData.irRaw{1}));
        for rep = 1:cfg.repetitions
            allRecorded = allRecorded + spkData.irRaw{rep};
        end
        avgRecorded = allRecorded / cfg.repetitions;
        
        % ============== 延迟估计 ==============
        % 采用三阶段估计策略
        
        % 阶段1: 基于物理约束和已知信息的初始估计
        fprintf('    阶段1: 基于物理约束的初始估计\n');
        preSilenceSamples = round(cfg.preSilenceSec * cfg.fs);
        deviceLatency = cfg.deviceLatencySamples; % 硬件延迟，需要预先标定
        
        % 计算理论最小延迟
        minTheoreticalDelay = preSilenceSamples + deviceLatency + cfg.minPhysDelaySamples;
        
        % 提取合适的信号段进行延迟估计
        % 从理论最小延迟前100ms开始，避免错过早期到达
        searchStart = max(1, minTheoreticalDelay - round(0.1 * cfg.fs));
        searchEnd = min(size(avgRecorded,1), searchStart + length(sweepDrive) * 2);
        
        recordedSegment = avgRecorded(searchStart:searchEnd, 1);
        
        % 阶段2: 多算法融合的精确延迟估计
        fprintf('    阶段2: 多算法融合精确估计\n');
        [delayEstimate_raw, delayMetrics] = estimate_delay_industrial_advanced(...
            recordedSegment, sweepDrive, cfg);
        
        % 阶段3: 修正为全局延迟并验证
        delayEstimate_global = delayEstimate_raw + (searchStart - 1);
        
        % 物理合理性检查
        [delayEstimate, validationResult] = validate_delay_estimate(...
            delayEstimate_global, avgRecorded(:,1), sweepDrive, cfg);
        
        fprintf('    延迟估计结果:\n');
        fprintf('      原始估计: %d 样本\n', delayEstimate_raw);
        fprintf('      全局延迟: %d 样本 (%.3f ms)\n', ...
            delayEstimate, delayEstimate/cfg.fs*1000);
        fprintf('      物理检查: %s\n', ternary(validationResult.pass, '通过', '失败'));
        fprintf('      置信度: %.1f%%\n', validationResult.confidence * 100);
        
        % 保存延迟估计的元数据
        spkData.delayEstimate = delayEstimate;
        spkData.delayMetrics = delayMetrics;
        spkData.delayValidation = validationResult;
        
        % === 阶段3: 精确反卷积 ===
        fprintf('[measure] 阶段3: 脉冲响应提取\n');
        
        for rep = 1:cfg.repetitions
            fprintf('  处理重复 %d... ', rep);
            
            % 确定截取窗口（增强鲁棒性）
            extractStart = max(1, delayEstimate - round(0.1*cfg.fs)); % 提前100ms
            extractEnd = min(size(spkData.irRaw{rep},1), ...
                extractStart + sweepCoreLen + Lh*2);
            
            if extractEnd - extractStart < sweepCoreLen + Lh
                extractStart = max(1, extractEnd - (sweepCoreLen + Lh*2));
            end
            
            % 对每个麦克风进行反卷积
            for m = 1:numErrMics
                % 截取信号段
                recSegment = spkData.irRaw{rep}(extractStart:extractEnd, m);
                
                % 反卷积（内部使用能量比SNR）
                irResult = deconv_industrial(recSegment, sweepDrive, cfg);
                
                % 保存结果
                irLength = min(Lh, length(irResult.ir));
                spkData.irFinal(1:irLength, m, rep) = irResult.ir(1:irLength);
                spkData.snrEst(rep, m) = irResult.snr;
                spkData.peakPos(rep, m) = irResult.peakPos;
            end
            
            fprintf('完成\n');
        end
        
        % === 阶段4: 数据质量评估 ===
        fprintf('[measure] 阶段4: 质量评估\n');
        
        % 添加相干性检查（对每个重复和麦克风）
        spkData.coherenceEst = zeros(cfg.repetitions, numErrMics);
        for rep = 1:cfg.repetitions
            for m = 1:numErrMics
                % 使用记录的原始信号和激励信号计算相干性
                rec_segment = spkData.irRaw{rep}(:, m);
                spkData.coherenceEst(rep, m) = compute_mean_coherence(...
                    rec_segment, sweepDrive, cfg.fs);
            end
        end

        % 计算统计指标
        qualityMetrics = assess_ir_quality(spkData.irFinal, spkData.snrEst, spkData.peakPos, spkData.coherenceEst, cfg);
        
        % 判断是否可用
        usable = qualityMetrics.usable;
        
        % === 阶段5: 平均与对齐 ===
        if usable
            fprintf('[measure] 阶段5: 数据平均\n');
            
            % 对齐所有重复
            irAligned = align_impulse_responses(spkData.irFinal, spkData.peakPos, cfg);
            
            % 加权平均（基于SNR）
            weights = mean(spkData.snrEst, 2);
            weights = weights - min(weights);
            if sum(weights) > 0
                weights = weights / sum(weights);
                irAvg = zeros(Lh, numErrMics);
                for rep = 1:cfg.repetitions
                    irAvg = irAvg + weights(rep) * squeeze(irAligned(:,:,rep));
                end
            else
                irAvg = mean(irAligned, 3);
            end
            
            % 最终IR处理（注意：不再归一化！）
            irAvg = postprocess_ir(irAvg, cfg);  % 此函数应移除归一化
            
            % 保存到主容器
            impulseResponses(:, :, spk) = irAvg;
            
            fprintf('    平均完成，SNR: %.1f dB, 稳定性: %.2f, 相干性: %.2f\n', ...
                qualityMetrics.medianSNR, qualityMetrics.stability, qualityMetrics.meanCoherence);
        else
            fprintf('[measure] 警告: 数据质量不足，使用零IR\n');
            fprintf('  原因: SNR=%.1f dB (阈值%.1f), 相干性=%.2f (阈值%.2f)\n', ...
                qualityMetrics.medianSNR, cfg.snrThresholdDB, ...
                qualityMetrics.meanCoherence, 0.7);
            impulseResponses(:, :, spk) = zeros(Lh, numErrMics);
        end
        
        % 保存元数据
        meta.perSpeaker{spk} = struct(...
            'delayEstimate', delayEstimate, ...
            'qualityMetrics', qualityMetrics, ...
            'usable', usable, ...
            'snrEst', spkData.snrEst, ...
            'peakPos', spkData.peakPos, ...
            'coherenceEst', spkData.coherenceEst, ... % 新增
            'extractWindow', [extractStart, extractEnd]);
        
        fprintf('[measure] 扬声器 %d 完成: %s\n', spk, ternary(usable, '可用', '不可用'));
    end
    
    %% ============== 构建输出结构 ==============
    fprintf('\n[output] 构建输出结构...\n');
    
    secondary = struct();
    secondary.impulseResponses = impulseResponses;
    secondary.fs = cfg.fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;
    secondary.errorMicPhysicalChannels = errMicIdx;
    secondary.irLength = Lh;
    secondary.description = 'Industrial-grade secondary path (v4.2 - optimized)';
    secondary.timestamp = datetime('now', 'TimeZone', 'local');
    secondary.measurementInfo = struct(...
        'sweepLength', totalSamples, ...
        'sweepCoreLength', sweepCoreLen, ...
        'repetitions', cfg.repetitions, ...
        'irMaxLen', Lh);
    
    % 计算推荐延迟
    delayVector = zeros(1, cfg.numSpeakers);
    for i = 1:cfg.numSpeakers
        if meta.perSpeaker{i}.usable
            delayVector(i) = meta.perSpeaker{i}.delayEstimate;
        else
            delayVector(i) = cfg.minPhysDelaySamples;
        end
    end
    secondary.delayEstimateSamples = delayVector;
    
    % 保存诊断信息
    if cfg.saveDiagnosticInfo
        secondary.meta = meta;
    end
    
    %% ============== 保存结果 ==============
    fprintf('[output] 保存结果...\n');
    
    % 创建目录
    saveDir = fileparts(cfg.secondaryPathFile);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
        fprintf('    创建目录: %s\n', saveDir);
    end
    
    % 保存数据
    save(cfg.secondaryPathFile, 'secondary', '-v7.3');
    fprintf('[output] 保存完成: %s\n', cfg.secondaryPathFile);
    
    %% ============== 生成质量报告 ==============
    if cfg.generateReport
        generate_quality_report(secondary, meta, cfg);
    end
    
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('[success] 次级路径测量完成!\n');
    fprintf('%s\n', repmat('=', 1, 60));
    
catch ME
    fprintf('\n[ERROR] 测量过程中发生错误:\n');
    fprintf('  消息: %s\n', ME.message);
    fprintf('  文件: %s\n', ME.stack(1).name);
    fprintf('  行号: %d\n', ME.stack(1).line);
    
    % 尝试保存已有数据
    if exist('secondary', 'var')
        try
            errorFile = strrep(cfg.secondaryPathFile, '.mat', '_ERROR.mat');
            save(errorFile, 'secondary', '-v7.3');
            fprintf('[ERROR] 部分数据已保存: %s\n', errorFile);
        catch
            fprintf('[ERROR] 无法保存数据\n');
        end
    elseif exist('meta', 'var')
        % 即使secondary不存在，也尝试保存meta数据
        try
            errorFile = strrep(cfg.secondaryPathFile, '.mat', '_META_ERROR.mat');
            save(errorFile, 'meta', '-v7.3');
            fprintf('[ERROR] 元数据已保存: %s\n', errorFile);
        catch
            fprintf('[ERROR] 无法保存元数据\n');
        end
    end
    
    rethrow(ME);
end

fprintf('\n[complete] 程序结束\n');

end % ← 主函数结束

%% ========================================================================
%% 子函数定义
%% ========================================================================

function cfg = validate_config(cfg)
% 配置验证与修正
requiredFields = {'fs', 'numSpeakers', 'irMaxLen', 'minPhysDelaySamples', ...
    'maxPhysDelaySamples', 'repetitions', 'snrThresholdDB'};

for i = 1:length(requiredFields)
    if ~isfield(cfg, requiredFields{i})
        error('配置缺少必要字段: %s', requiredFields{i});
    end
end

% 参数合理性检查
if cfg.minPhysDelaySamples < 1
    cfg.minPhysDelaySamples = 1;
    fprintf('[config] 警告: minPhysDelaySamples调整到1\n');
end

if cfg.maxPhysDelaySamples < cfg.minPhysDelaySamples + 100
    cfg.maxPhysDelaySamples = cfg.minPhysDelaySamples + 10000;
    fprintf('[config] 警告: maxPhysDelaySamples调整到%d\n', cfg.maxPhysDelaySamples);
end

if cfg.repetitions < 1
    cfg.repetitions = 1;
    fprintf('[config] 警告: repetitions调整到1\n');
end

% 设置默认值
if ~isfield(cfg, 'deconvRegEps')
    cfg.deconvRegEps = 1e-12;
end
if ~isfield(cfg, 'enableLowFreqBoost')
    cfg.enableLowFreqBoost = false;
end
if ~isfield(cfg, 'saveDiagnosticInfo')
    cfg.saveDiagnosticInfo = true;
end
if ~isfield(cfg, 'generateReport')
    cfg.generateReport = false;
end
if ~isfield(cfg, 'coherenceThreshold')
    cfg.coherenceThreshold = 0.7; % 新增默认相干性阈值
end
end

function str = ternary(condition, trueStr, falseStr)
% 三元操作符模拟
if condition
    str = trueStr;
else
    str = falseStr;
end
end