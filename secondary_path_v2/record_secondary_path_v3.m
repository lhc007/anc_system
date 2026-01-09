function record_secondary_path_v3()
% record_secondary_path_v3.m 多通道管道ANC次级路径测量系统
% 版本：v4.7（完全修复版）

clc;

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
% 全局硬件句柄（用于确保即使出错也能释放）
global_hw = [];

try
    %% ============== 生成激励信号 ==============
    fprintf('[signal] 生成ESS激励信号...\n');
    [sweepCore_scaled, sweepSig, exciteInfo, weightedSweepForCoherence, invFilter] = generate_excitation_signal(cfg);
    sweepCoreLen = exciteInfo.coreLength;
    fs = cfg.fs;    
    % ========== 信号验证 ==========
    totalSamples = length(sweepSig); 
    essTotalLen = exciteInfo.totalLength;
    
    fprintf('[signal] 信号参数:\n');
    fprintf('  - 核心扫频时长: %.3f s (%d 样本)\n', exciteInfo.T, sweepCoreLen);
    fprintf('  - ESS激励总长: %d 样本 (含 %d ms 前静音, %d ms 后静音)\n', ...
        essTotalLen, round(exciteInfo.padLeadingSamples*1000/fs), round(exciteInfo.padTrailingSamples*1000/fs));
    fprintf('  - 采样率: %d Hz\n', fs);
    
    if length(sweepCore_scaled) ~= sweepCoreLen
        error('扫频核心长度不匹配！');
    end
    
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
    
    %% ============== 生成静音信号用于硬件重置 ==============
    silence_duration = 0.5;  % 100ms静音用于重置缓冲区
    silence_signal = zeros(round(silence_duration * fs), 1);
    
    %% ============== 主测量循环 ==============
    fprintf('\n[measure] 开始次级路径测量...\n');
    hw = hardware_init_measure(cfg);
    global_hw = hw;  % 保存到全局变量用于错误处理

    for spk = 1:cfg.numSpeakers
        fprintf('\n%s\n', repmat('=', 1, 60));
        fprintf('[measure] 扬声器 %d/%d\n', spk, cfg.numSpeakers);
        fprintf('%s\n', repmat('=', 1, 60));
                
        % ========== 准备激励信号 ==========
        spkDriveSig = cfg.spkAmplitude(spk) * sweepSig(:);
        maxVal = max(abs(spkDriveSig));
        if maxVal > 0.95
            scale = 0.9 / maxVal;
            spkDriveSig = spkDriveSig * scale;
        else
            % 确保最小幅度（对抗环境噪声）
            minAmp = 0.3;
            if maxVal < minAmp
                spkDriveSig = spkDriveSig * (minAmp / maxVal);
                fprintf('  提升激励幅度至 %.2f\n', minAmp);
            end
        end
        % ========== 预播放静音，稳定硬件 ==========
        fprintf('  预播放静音稳定硬件...\n');
        [recordedFull, info] = play_and_record(hw, spkDriveSig, spk, cfg);
        save('debug_recorded.mat', 'recordedFull', 'spkDriveSig', 'cfg'); % 保存原始数据
        fprintf('  播放长度: %d, 录制长度: %d\n', length(spkDriveSig), size(recordedFull,1));
       
        % ========== 数据采集容器 ==========
        irFinal = zeros(Lh, numErrMics, cfg.repetitions);
        snrEst = zeros(cfg.repetitions, numErrMics);
        peakPos = zeros(cfg.repetitions, numErrMics);
        coherenceEst = zeros(cfg.repetitions, numErrMics);
        recordedRaw = cell(cfg.repetitions, 1);  % 存储原始录制数据用于调试
        
        fprintf('[measure] 阶段1: 数据采集 (%d次重复)\n', cfg.repetitions);
        
        for rep = 1:cfg.repetitions
            fprintf('  重复 %d/%d... ', rep, cfg.repetitions);
            
            % ========== 播放前静音重置缓冲区 ==========
            if rep > 1
                [~, ~] = play_and_record(hw, silence_signal, spk, cfg);
            end
            
            % ========== 播放并录制主信号 ==========
            [recordedFull, ~] = play_and_record(hw, spkDriveSig, spk, cfg);
            recordedRaw{rep} = recordedFull;  % 保存原始数据
            
            % ========== 对齐信号 ==========
            [actual_start_idx, recorded_aligned, clickAbsPos] = align_sweep_start(recordedFull, exciteInfo, errMicIdx);
           
            % figure;
            % subplot(2,1,1);
            % plot(recorded_aligned(:,1));
            % title('对齐后的录制信号（通道1）');
            % xlabel('样本'); ylabel('幅度');
            % 
            % subplot(2,1,2);
            % N = length(recorded_aligned(:,1));
            % f = (0:N-1)*(fs/N);
            % Y = fft(recorded_aligned(:,1));
            % semilogx(f(1:N/2), 20*log10(abs(Y(1:N/2))));
            % title('频谱（应看到 100–4000 Hz 扫频能量）');
            % xlabel('频率 (Hz)'); ylabel('幅度 (dB)');
            % grid on;
            % ========== 反卷积处理 ==========
            for m = 1:numErrMics
                recSegment = recorded_aligned(:, m);
                
                % 反卷积
                irResult = deconv_ess(recSegment, invFilter, cfg, clickAbsPos, exciteInfo.coreLength);

                 % 截取到指定长度
                irLength = min(Lh, length(irResult.ir));
                irFinal(1:irLength, m, rep) = irResult.ir(1:irLength);
                snrEst(rep, m) = irResult.snr;
                peakPos(rep, m) = irResult.peakInSegment;  % 应为常数（如 128）
                
                % 计算相干性     
                coherenceEst(rep, m) = compute_reconstruction_coherence(...
                    irResult.ir, ...
                    weightedSweepForCoherence, ...
                    recSegment);
            end
            fprintf('完成 (SNR: %.1f dB, 延迟: %d)\n', ...
                mean(snrEst(rep, :)),  round(median(peakPos(rep, :))));
            
            % ========== 重复间隔 ==========
            if rep < cfg.repetitions
                pause(cfg.repetitionInterval);
            end
        end
         % ========== 硬件释放 ==========
        if ~isempty(hw)
            hw.release();
            global_hw = [];
        end

        % ========== 数据质量评估 ==========
        fprintf('[measure] 阶段2: 质量评估\n');
        metrics = assess_ir_quality(irFinal, snrEst, peakPos, coherenceEst, cfg);

        usable = metrics.usable;
        
        % ========== 延迟估计（统一策略）==========
        % 新语法：仅支持MATLAB R2018b及以上版本 新语法更灵活，可以指定维度
        medianPeakPos = median(peakPos, 'all');  

        % 物理延迟 = 峰值位置 + 硬件延迟修正 - 1
        % （峰值位置是从1开始的索引，实际延迟是峰值位置-1）
        if isfield(cfg, 'hardwareDelaySamples')
            hardwareDelay = cfg.hardwareDelaySamples;
        else
            hardwareDelay = 0;  % 默认为0，需要校准
        end
        
        delayEstimate = (medianPeakPos - 1) + hardwareDelay;
        
        % 确保延迟在合理范围内
        delayEstimate = max(cfg.minPhysDelaySamples, ...
            min(cfg.maxPhysDelaySamples, round(delayEstimate)));
        
        %% =================================
        % 校准值（通过 loopback 实验测得）
         hardwareDelayCalibrated = 120;  % 单位：samples
         
        % 计算物理延迟（去除系统固有延迟）
        delaySamples = (metrics.medianPeakPosPerMic - 1) - hardwareDelayCalibrated;
         delaySeconds = delaySamples / cfg.fs;
        
        fprintf('各通道物理延迟（秒）: %s\n', mat2str(delaySeconds'));
        %% =================================

        % ========== 数据处理 ==========
        if usable
            fprintf('[measure] 阶段3: 数据处理\n');
            
            % 对齐IR
            irAligned = align_impulse_responses(irFinal, peakPos, cfg);
            
            % SNR加权平均
            weights = mean(snrEst, 2);
            weights = weights - min(weights);
            
            if sum(weights) > 0 && cfg.repetitions > 1
                weights = weights / sum(weights);
                irAvg = zeros(Lh, numErrMics);
                for rep = 1:cfg.repetitions
                    irAvg = irAvg + weights(rep) * squeeze(irAligned(:,:,rep));
                end
            else
                irAvg = mean(irAligned, 3);
            end
            
            % 后处理
            irAvg = postprocess_ir(irAvg, cfg);
            impulseResponses(:, :, spk) = irAvg;
            
            fprintf('  处理完成: SNR=%.1fdB, 相干性=%.2f, 延迟=%d样本\n', ...
                metrics.medianSNR, metrics.meanCoherence, delayEstimate);
        else
            fprintf('[measure] 警告: 数据质量不足，使用零IR\n');
            impulseResponses(:, :, spk) = zeros(Lh, numErrMics);
        end
        
        % ========== 保存元数据 ==========
        % 创建标量结构体
        speakerMeta = struct();
        
        speakerMeta.delayEstimate = delayEstimate;
        speakerMeta.medianPeakPos = medianPeakPos;
        speakerMeta.hardwareDelay = hardwareDelay;
        speakerMeta.qualityMetrics = metrics;
        speakerMeta.usable = usable;
        speakerMeta.snrEst = snrEst;
        speakerMeta.peakPos = peakPos;
        speakerMeta.coherenceEst = coherenceEst;
        speakerMeta.rawRecordings = recordedRaw;  % cell array 是合法字段值！
        
        meta.perSpeaker{spk} = speakerMeta;  % 明确赋值为标量结构体

        fprintf('[measure] 扬声器 %d 完成: %s\n', ...
            spk, ternary(usable, '✅ 可用', '❌ 不可用'));
        
        % ========== 扬声器间暂停（避免串扰）==========
        if spk < cfg.numSpeakers
            fprintf('  扬声器间暂停 %.1f 秒...\n', cfg.repetitionInterval * 2);
            pause(cfg.repetitionInterval * 2);
        end
    end
    
    %% ============== 构建输出结构 ==============
    secondary = struct();
    secondary.impulseResponses = impulseResponses;
    secondary.fs = cfg.fs;
    secondary.numSpeakers = cfg.numSpeakers;
    secondary.numMics = numErrMics;
    secondary.errorMicPhysicalChannels = errMicIdx;
    secondary.irLength = Lh;
    secondary.description = 'Industrial-grade secondary path (v4.7)';
    secondary.timestamp = datetime('now', 'TimeZone', 'local');
    
    % 延迟向量
    delayVector = zeros(1, cfg.numSpeakers);
    for i = 1:cfg.numSpeakers
        if meta.perSpeaker{i}.usable
            delayVector(i) = meta.perSpeaker{i}.delayEstimate;
        else
            delayVector(i) = cfg.minPhysDelaySamples;
        end
    end
    secondary.delayEstimateSamples = delayVector;
    
    % 可选：保存诊断信息
    if cfg.saveDiagnosticInfo
        secondary.meta = meta;
        % 添加配置摘要（不保存可能包含敏感信息的完整配置）
        secondary.configSummary = struct(...
            'fs', cfg.fs, ...
            'repetitions', cfg.repetitions, ...
            'irMaxLen', cfg.irMaxLen, ...
            'numSpeakers', cfg.numSpeakers, ...
            'numErrorMics', numErrMics);
    end
    
    %% ============== 保存结果 ==============
    saveDir = fileparts(cfg.secondaryPathFile);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    % 添加版本信息
    secondary.version = 'v4.7';
    secondary.checksum = compute_data_checksum(impulseResponses);
    
    save(cfg.secondaryPathFile, 'secondary', '-v7.3');
    fprintf('[output] 保存完成: %s\n', cfg.secondaryPathFile);
    
    if cfg.generateReport
        generate_quality_report(secondary, meta, cfg);
    end
    
    fprintf('\n[success] 次级路径测量完成!\n');
    
catch ME
    % ========== 错误处理：确保硬件释放 ==========
    if ~isempty(global_hw)
       global_hw.release();
       
    end
    fprintf('\n[ERROR] 测量失败: %s\n', ME.message);
    fprintf('堆栈跟踪:\n');
    for i = 1:min(5, length(ME.stack))  % 只显示前5个堆栈帧
        fprintf('  %s (第%d行)\n', ME.stack(i).name, ME.stack(i).line);
    end
    
    % 保存部分数据用于调试
    try
        errorData = struct(...
            'errorMessage', ME.message, ...
            'errorTime', datetime('now', 'TimeZone', 'local'), ...
            'partialData', impulseResponses, ...
            'lastSpeaker', spk);
        
        errorFile = sprintf('secondary_path_error_%s.mat', ...
            datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        save(errorFile, 'errorData', '-v7.3');
        fprintf('[output] 错误数据已保存: %s\n', errorFile);
    catch
        fprintf('  无法保存错误数据\n');
    end
    
    rethrow(ME);
end

fprintf('\n[complete] 程序结束\n');
end

%% ========================================================================
% 修复的辅助函数
% ========================================================================

function cfg = validate_config(cfg)
requiredFields = {'fs', 'numSpeakers', 'irMaxLen', 'minPhysDelaySamples', ...
    'maxPhysDelaySamples', 'repetitions', 'snrThresholdDB', 'micChannels'};
for i = 1:length(requiredFields)
    if ~isfield(cfg, requiredFields{i})
        error('配置缺少必要字段: %s', requiredFields{i});
    end
end

% 自动补充缺失字段
if ~isfield(cfg, 'hardwareDelaySamples')
    cfg.hardwareDelaySamples = 0;  % 需要实际测量校准
end

if ~isfield(cfg, 'repetitionInterval')
    cfg.repetitionInterval = 0.5;
end

if ~isfield(cfg, 'saveDiagnosticInfo')
    cfg.saveDiagnosticInfo = false;
end

if ~isfield(cfg, 'generateReport')
    cfg.generateReport = false;
end

% 验证值范围
if cfg.minPhysDelaySamples < 1
    cfg.minPhysDelaySamples = 1;
end

if cfg.maxPhysDelaySamples < cfg.minPhysDelaySamples + 50
    cfg.maxPhysDelaySamples = cfg.minPhysDelaySamples + 200;
end

if cfg.repetitions < 1
    cfg.repetitions = 1;
end

if ~isfield(cfg, 'deconvPreDelayKeep')
    cfg.deconvPreDelayKeep = 64;
end

% 验证麦克风通道
if ~isfield(cfg.micChannels, 'error')
    error('配置缺少误差麦克风通道设置');
end
end

function checksum = compute_data_checksum(data)
% 计算数据校验和，用于验证数据完整性
data_flat = data(:);
checksum = sum(abs(data_flat(1:min(1000, length(data_flat)))));
checksum = mod(checksum, 2^31);  % 避免溢出
end