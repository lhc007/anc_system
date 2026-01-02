function record_secondary_path_v3()
% record_secondary_path_v3.m 多通道管道ANC次级路径测量系统
% 版本：v4.5

clear; clc;

%% ============== 配置加载与验证 ==============
cfg = anc_config();
fs = cfg.fs();

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
    [sweepCore_scaled, sweepSig, info] = generate_excitation_signal(cfg);
    
    % ========== 信号验证 ==========
    sweepCoreLen = info.N;                  % 核心扫频长度
    essTotalLen = info.totalLength;         % ESS总长度
    totalSamples = length(sweepSig); 
    
    fprintf('[signal] 信号参数:\n');
    fprintf('  - 核心扫频时长: %.3f s (%d 样本)\n', info.sweepDuration, sweepCoreLen);
    fprintf('  - ESS激励总长: %d 样本 (含 %d ms 前静音, %d ms 后静音)\n', ...
        essTotalLen, round(info.padLeading*1000), round(info.padTrailing*1000));
    fprintf('  - 采样率: %d Hz\n', fs);
    
    % 验证扫频信号长度
    if length(sweepCore_scaled) ~= sweepCoreLen
        error('扫频核心长度不匹配！sweepCore_scaled长度为%d，应为%d', ...
            length(sweepCore_scaled), sweepCoreLen);
    end
    fprintf('  [验证] sweepCore_scaled长度: %d (正确)\n', length(sweepCore_scaled));

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
        spkData.irAligned = cell(cfg.repetitions, 1);
        spkData.delayEst = zeros(cfg.repetitions, 1);
        spkData.snrEst = zeros(cfg.repetitions, numErrMics);
        spkData.peakPos = zeros(cfg.repetitions, numErrMics);
        spkData.coherenceEst = zeros(cfg.repetitions, numErrMics);
        spkData.irFinal = zeros(Lh, numErrMics, cfg.repetitions);
        
        % === 阶段1: 数据采集（基于扫频起始点对齐） ===
        fprintf('[measure] 阶段1: 数据采集（基于扫频起始点对齐）\n');
        
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
            [recordedFull, ~] = play_and_record(hw, spkDriveSig, spk, cfg);
            
            % 基于扫频起始点对齐
            [actual_start_idx, recorded_aligned] = align_sweep_start(...
                recordedFull, errMicIdx, sweepCore_scaled, cfg);
            
            % 保存对齐后数据
            spkData.irAligned{rep} = recorded_aligned;
            
            % 释放硬件
            hw.release();
            
            fprintf('完成 (扫频起始@%d, 对齐信号长度: %d)\n', ...
                actual_start_idx, size(recorded_aligned,1));
        end

        %% === 阶段2: 自适应延迟估计（基于对齐信号，多通道） ===
        fprintf('[measure] 阶段2: 自适应延迟估计（基于对齐信号，多通道）\n');
        
        % 使用所有误差麦克风通道进行延迟估计
        [allDelays, allConfidences, allCorrelations] = estimate_delay_multichannel(...
            spkData.irAligned, sweepCore_scaled, cfg, numErrMics);
        
        % 使用高置信度通道的延迟（或加权平均）
        highConfMask = allConfidences > 0.5; % 置信度阈值
        if sum(highConfMask) >= 1
            % 使用高置信度通道的加权平均
            weights = allConfidences(highConfMask);
            if sum(weights) > 0
                weights = weights / sum(weights);
                delayEstimate = round(sum(allDelays(highConfMask) .* weights));
            else
                delayEstimate = round(mean(allDelays));
            end
            fprintf('  使用加权平均延迟: %d 样本\n', delayEstimate);
        else
            % 所有通道置信度都低，使用第一个通道
            delayEstimate = allDelays(1);
            fprintf('  警告: 所有通道置信度低，使用通道1延迟: %d 样本\n', delayEstimate);
        end
        
        % 验证结果（基于最终选择的延迟）
        delayStd = std(allDelays);
        delayMean = mean(allDelays);
        delayRange = max(allDelays) - min(allDelays);
        fprintf('  延迟统计: 平均=%.1f, 标准差=%.1f, 范围=[%d,%d]\n', ...
            delayMean, delayStd, min(allDelays), max(allDelays));
        
        % 更严格的延迟范围检查（5ms）
        if delayRange > round(0.005 * fs) % >5ms差异
            fprintf('  警告: 通道间延迟差异较大 (%.1f ms)，可能影响精度\n', delayRange/fs*1000);
            validationResult = struct('pass', false, 'confidence', 0.6, ...
                'delay_global', delayEstimate, 'correlation', mean(allCorrelations));
        else
            validationResult = struct('pass', true, 'confidence', min(0.9, mean(allConfidences)), ...
                'delay_global', delayEstimate, 'correlation', mean(allCorrelations));
        end
        
        % 保存延迟估计的元数据
        spkData.delayEstimate = delayEstimate;
        spkData.delayValidation = validationResult;
        spkData.delayPerChannel = allDelays;
        spkData.delayConfidence = allConfidences;
        spkData.delayCorrelation = allCorrelations;
        
        % === 阶段3: 精确反卷积（使用对齐信号） ===
        fprintf('[measure] 阶段3: 脉冲响应提取（使用对齐信号）\n');
        
        for rep = 1:cfg.repetitions
            fprintf('  处理重复 %d... ', rep);
            
            % 使用对齐后的信号进行反卷积
            for m = 1:numErrMics
                % 对齐后的信号段
                recSegment = spkData.irAligned{rep}(:, m);
                fprintf('[DEBUG] 反卷积参数: irMaxLen=%d, regEps=%.1e\n', ...
                    cfg.irMaxLen, cfg.deconvRegEps);
                % 反卷积
                irResult = deconv_industrial(recSegment, sweepCore_scaled, cfg);
                
                % 保存IR结果
                irLength = min(Lh, length(irResult.ir));
                spkData.irFinal(1:irLength, m, rep) = irResult.ir(1:irLength);
                spkData.snrEst(rep, m) = irResult.snr;
                spkData.peakPos(rep, m) = irResult.peakPos;
                
                % 修复：使用重建相干性计算方法
                % 比较实际录制信号与通过提取IR重建的信号
                spkData.coherenceEst(rep, m) = compute_reconstruction_coherence(...
                    recSegment, sweepCore_scaled, irResult.ir, fs);
            end
            
            fprintf('完成\n');
        end

        % 计算统计指标
        qualityMetrics = assess_ir_quality(spkData.irFinal, spkData.snrEst, ...
            spkData.peakPos, spkData.coherenceEst, cfg);
        
        % 判断是否可用
        usable = qualityMetrics.usable;
        
        % === 阶段5: 数据平均 ===
        if usable
            fprintf('[measure] 阶段5: 数据平均\n');
            
            % 对齐所有重复（基于峰值位置）
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
            
            % 最终IR处理
            irAvg = postprocess_ir(irAvg, cfg);
            
            % 保存到主容器
            impulseResponses(:, :, spk) = irAvg;
            
            fprintf('    平均完成，SNR: %.1f dB, 稳定性: %.2f, 相干性: %.2f\n', ...
                qualityMetrics.medianSNR, qualityMetrics.stability, qualityMetrics.meanCoherence);
        else
            fprintf('[measure] 警告: 数据质量不足，使用零IR\n');
            fprintf('  原因:\n');
            if ~qualityMetrics.snrOK
                fprintf('    - SNR=%.1f dB < 阈值%.1f dB\n', ...
                    qualityMetrics.medianSNR, cfg.snrThresholdDB);
            end
            if ~qualityMetrics.coherenceOK
                fprintf('    - 相干性=%.2f < 阈值%.2f\n', ...
                    qualityMetrics.medianCoherence, cfg.coherenceThreshold);
            end
            if ~qualityMetrics.stablePeaks && cfg.repetitions > 1
                fprintf('    - 峰值标准差=%.1f > 阈值%.1f\n', ...
                    qualityMetrics.peakStd, cfg.maxPeakStd);
            end
            impulseResponses(:, :, spk) = zeros(Lh, numErrMics);
        end
        
        % 保存元数据
        meta.perSpeaker{spk} = struct(...
            'delayEstimate', delayEstimate, ...
            'delayPerChannel', spkData.delayPerChannel, ...
            'delayConfidence', spkData.delayConfidence, ...
            'qualityMetrics', qualityMetrics, ...
            'usable', usable, ...
            'snrEst', spkData.snrEst, ...
            'peakPos', spkData.peakPos, ...
            'coherenceEst', spkData.coherenceEst, ...
            'alignmentMethod', 'sweep_start_point', ...
            'coherenceMethod', 'reconstruction_coherence'); % 新增：记录相干性计算方法
        
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
    secondary.description = 'Industrial-grade secondary path (v4.5 - 修复相干性计算)';
    secondary.timestamp = datetime('now', 'TimeZone', 'local');
    secondary.measurementInfo = struct(...
        'sweepLength', length(sweepSig), ...
        'sweepCoreLength', sweepCoreLen, ...
        'repetitions', cfg.repetitions, ...
        'irMaxLen', Lh, ...
        'coherenceMethod', 'reconstruction_based'); % 新增：记录相干性计算方法
    
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
    
    % 打印最终统计
    fprintf('\n[统计] 测量结果汇总:\n');
    usableCount = 0;
    for i = 1:cfg.numSpeakers
        if meta.perSpeaker{i}.usable
            usableCount = usableCount + 1;
            fprintf('  扬声器 %d: SNR=%.1f dB, 相干性=%.2f, 延迟=%d样本\n', ...
                i, meta.perSpeaker{i}.qualityMetrics.medianSNR, ...
                meta.perSpeaker{i}.qualityMetrics.medianCoherence, ...
                meta.perSpeaker{i}.delayEstimate);
        else
            fprintf('  扬声器 %d: 不可用\n', i);
        end
    end
    fprintf('  可用率: %d/%d (%.1f%%)\n', usableCount, cfg.numSpeakers, ...
        usableCount/cfg.numSpeakers*100);
    
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

end 

%% ========================================================================
% 子函数定义
% ========================================================================

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
    cfg.maxPhysDelaySamples = cfg.minPhysDelaySamples + 100000;
    fprintf('[config] 警告: maxPhysDelaySamples调整到%d\n', cfg.maxPhysDelaySamples);
end

if cfg.repetitions < 1
    cfg.repetitions = 1;
    fprintf('[config] 警告: repetitions调整到1\n');
end

if ~isfield(cfg, 'deconvPreDelayKeep')
    cfg.deconvPreDelayKeep = 10000;
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