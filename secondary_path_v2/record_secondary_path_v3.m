function record_secondary_path_v3()
% record_secondary_path_v3.m 多通道管道ANC次级路径测量系统
% 版本：v4.5（修复版）
% sweepSig：完整 ESS 信号（含前后静音），用于播放和对齐。
% sweepCore_scaled：仅核心扫频部分（无静音），用于反卷积。

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
    sweepCoreLen = info.N;                  
    essTotalLen = info.totalLength;         
    totalSamples = length(sweepSig); 
    
    fprintf('[signal] 信号参数:\n');
    fprintf('  - 核心扫频时长: %.3f s (%d 样本)\n', info.sweepDuration, sweepCoreLen);
    fprintf('  - ESS激励总长: %d 样本 (含 %d ms 前静音, %d ms 后静音)\n', ...
        essTotalLen, round(info.padLeading*1000), round(info.padTrailing*1000));
    fprintf('  - 采样率: %d Hz\n', fs);
    
    if length(sweepCore_scaled) ~= sweepCoreLen
        error('扫频核心长度不匹配！');
    end
    fprintf('  [验证] sweepCore_scaled长度: %d (正确)\n', length(sweepCore_scaled));

    %% ============== 容器初始化 ==============
    Lh = cfg.irMaxLen;
    impulseResponses = zeros(Lh, numErrMics, cfg.numSpeakers);
    
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
        
        spkData = struct();
        spkData.irRaw = cell(cfg.repetitions, 1);
        spkData.irAligned = cell(cfg.repetitions, 1);
        spkData.snrEst = zeros(cfg.repetitions, numErrMics);
        spkData.peakPos = zeros(cfg.repetitions, numErrMics);
        spkData.coherenceEst = zeros(cfg.repetitions, numErrMics);
        spkData.irFinal = zeros(Lh, numErrMics, cfg.repetitions);
        
        % === 阶段1: 数据采集 ===
        fprintf('[measure] 阶段1: 数据采集\n');
        
        for rep = 1:cfg.repetitions
            fprintf('  重复 %d/%d... ', rep, cfg.repetitions);
            
            if rep == 1
                hw = hardware_init_measure(cfg);
                fprintf('硬件已初始化\n');
            else
                pause(0.5);
                hw = hardware_init_measure(cfg);
            end
            
            spkDriveSig = cfg.spkAmplitude(spk) * sweepSig(:);
            maxVal = max(abs(spkDriveSig));
            if maxVal > 0.95
                scale = 0.9 / maxVal;
                spkDriveSig = spkDriveSig * scale;
                fprintf('    防削波缩放: %.2f\n', scale);
            end
            
            % [播放 & 录制]
            [recordedFull, ~] = play_and_record(hw, spkDriveSig, spk, cfg);
            
            % ✅ 修复1: 使用 recordedFull，不是 recorded_raw！
            [~, recorded_aligned] = align_sweep_start(recordedFull, sweepSig, cfg, errMicIdx);
            
            spkData.irAligned{rep} = recorded_aligned;
            hw.release();
            fprintf('完成\n');
        end

        % === 阶段2: 反卷积 + 相干性计算（跳过冗余延迟估计）===
        fprintf('[measure] 阶段2: 冲激响应提取与质量评估\n');
        
        for rep = 1:cfg.repetitions
            fprintf('  处理重复 %d... ', rep);
            for m = 1:numErrMics
                recSegment = spkData.irAligned{rep}(:, m);
                
                % 反卷积
                irResult = deconv_industrial(recSegment, sweepCore_scaled, cfg);
                
                irLength = min(Lh, length(irResult.ir));
                spkData.irFinal(1:irLength, m, rep) = irResult.ir(1:irLength);
                spkData.snrEst(rep, m) = irResult.snr;
                spkData.peakPos(rep, m) = irResult.peakPos;
                
                % ✅ 修复2: 正确调用 compute_reconstruction_coherence
                % 参数顺序: (recorded, ir, sweep_ref, peakPos, fs)
                spkData.coherenceEst(rep, m) = compute_reconstruction_coherence(...
                    recSegment, irResult.ir, sweepCore_scaled, irResult.peakPos, fs);
            end
            fprintf('完成\n');
        end

        % 质量评估
        qualityMetrics = assess_ir_quality(spkData.irFinal, spkData.snrEst, ...
            spkData.peakPos, spkData.coherenceEst, cfg);
        usable = qualityMetrics.usable;
        
        % === 阶段3: 数据平均 ===
        if usable
            fprintf('[measure] 阶段3: 数据平均\n');
            irAligned = align_impulse_responses(spkData.irFinal, spkData.peakPos, cfg);
            
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
            
            irAvg = postprocess_ir(irAvg, cfg);
            impulseResponses(:, :, spk) = irAvg;
            
            fprintf('    平均完成，SNR: %.1f dB, 相干性: %.2f\n', ...
                qualityMetrics.medianSNR, qualityMetrics.meanCoherence);
        else
            fprintf('[measure] 警告: 数据质量不足，使用零IR\n');
            impulseResponses(:, :, spk) = zeros(Lh, numErrMics);
        end
        
        % ✅ 修复3: 使用 peakPos 计算真实物理延迟！
        avgPeakPos = mean(spkData.peakPos(:));
        delayEstimate = avgPeakPos - 1;  % ←←← 黄金标准！
        
        meta.perSpeaker{spk} = struct(...
            'delayEstimate', delayEstimate, ...   % ← 正确延迟
            'qualityMetrics', qualityMetrics, ...
            'usable', usable, ...
            'snrEst', spkData.snrEst, ...
            'peakPos', spkData.peakPos, ...
            'coherenceEst', spkData.coherenceEst);
        
        fprintf('[measure] 扬声器 %d 完成: %s (延迟=%d样本)\n', ...
            spk, ternary(usable, '可用', '不可用'), delayEstimate);
    end
    
    %% ============== 构建输出结构 ==============
    secondary = struct();
    secondary.impulseResponses = impulseResponses;
    secondary.fs = cfg.fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;
    secondary.errorMicPhysicalChannels = errMicIdx;
    secondary.irLength = Lh;
    secondary.description = 'Industrial-grade secondary path (v4.5 - 修复版)';
    secondary.timestamp = datetime('now', 'TimeZone', 'local');
    
    % ✅ 关键：使用基于 peakPos 的延迟
    delayVector = zeros(1, cfg.numSpeakers);
    for i = 1:cfg.numSpeakers
        if meta.perSpeaker{i}.usable
            delayVector(i) = meta.perSpeaker{i}.delayEstimate;  % ← 正确值
        else
            delayVector(i) = cfg.minPhysDelaySamples;
        end
    end
    secondary.delayEstimateSamples = delayVector;
    
    if cfg.saveDiagnosticInfo
        secondary.meta = meta;
    end
    
    %% ============== 保存结果 ==============
    saveDir = fileparts(cfg.secondaryPathFile);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    save(cfg.secondaryPathFile, 'secondary', '-v7.3');
    fprintf('[output] 保存完成: %s\n', cfg.secondaryPathFile);
    
    if cfg.generateReport
        generate_quality_report(secondary, meta, cfg);
    end
    
    fprintf('\n[success] 次级路径测量完成!\n');
    
catch ME
    fprintf('\n[ERROR] %s\n', ME.message);
    rethrow(ME);
end

fprintf('\n[complete] 程序结束\n');
end 

%% ========================================================================
% 子函数（保持不变）
% ========================================================================

function cfg = validate_config(cfg)
requiredFields = {'fs', 'numSpeakers', 'irMaxLen', 'minPhysDelaySamples', ...
    'maxPhysDelaySamples', 'repetitions', 'snrThresholdDB'};
for i = 1:length(requiredFields)
    if ~isfield(cfg, requiredFields{i})
        error('配置缺少必要字段: %s', requiredFields{i});
    end
end
if cfg.minPhysDelaySamples < 1
    cfg.minPhysDelaySamples = 1;
end
if cfg.maxPhysDelaySamples < cfg.minPhysDelaySamples + 100
    cfg.maxPhysDelaySamples = cfg.minPhysDelaySamples + 1000;
end
if cfg.repetitions < 1
    cfg.repetitions = 1;
end
if ~isfield(cfg, 'deconvPreDelayKeep')
    cfg.deconvPreDelayKeep = 64; % 默认保留64样本前置
end
end

function str = ternary(condition, trueStr, falseStr)
if condition
    str = trueStr;
else
    str = falseStr;
end
end