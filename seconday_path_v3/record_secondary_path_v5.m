function record_secondary_path_v5()
% record_secondary_path_v5.m —— 多通道ANC次级路径测量系统（v5.3 - Unified Coherence & Logging）
% 目标：为FXLMS提供高保真、LTI验证的次级路径模型
% 改进：
%   - 修复 coherent 计算接口不一致问题
%   - 统一变量命名（y_rec_crop）
%   - 增强日志输出（激励构成、对齐细节）
%   - 添加配置字段默认值
clc;
fprintf('[init] 启动次级路径测量系统 (v5.3 - Unified Coherence & Logging)...\n');

%% ============== 配置加载与验证 ==============
try
    cfg = anc_config();
catch
    error('无法加载 anc_config.m 文件');
end

cfg = validate_config_final(cfg);

errMicIdx = cfg.micChannels.error;
numErrMics = length(errMicIdx);
fs = cfg.fs;
Lh = cfg.irMaxLen;

fprintf('[init] 采样率: %d Hz | 扬声器数: %d | 误差麦通道: %s\n', ...
    fs, cfg.numSpeakers, mat2str(errMicIdx));

%% ============== 生成激励信号 ==============
[sweepSig, exciteInfo] = generate_excitation_signal(cfg);
% === 完善激励信号日志 ===
fprintf('[signal] 激励信号构成:\n');
fprintf('  前导静音: %.2f s (%d 样本)\n', exciteInfo.preSilence / fs, exciteInfo.preSilence);
fprintf('  核心扫频: %.2f s (%d 样本)\n', exciteInfo.T, exciteInfo.coreLength);
fprintf('  后导静音: %.2f s (%d 样本)\n', exciteInfo.postSilence / fs, exciteInfo.postSilence);
fprintf('  总长度  : %.2f s (%d 样本)\n', exciteInfo.totalLength / fs, exciteInfo.totalLength);

%% ============== 初始化 ==============
impulseResponses = zeros(Lh, numErrMics, cfg.numSpeakers);
global_hw = [];
allSNRs_all = zeros(cfg.repetitions, numErrMics, cfg.numSpeakers);
allCohs_all = zeros(cfg.repetitions, numErrMics, cfg.numSpeakers);

try
    hw = hardware_init_measure(cfg);
    global_hw = hw;

    %% ============== 主测量循环 ==============
    for spk = 1:cfg.numSpeakers
        fprintf('\n[measure] ▶ 扬声器 %d/%d\n', spk, cfg.numSpeakers);
        
        % 准备激励信号
        amp = cfg.spkAmplitude(spk);
        spkSig = amp * sweepSig(:);
        spkSig = safe_normalize(spkSig, cfg);

        % 预静音清缓冲
        long_silence = zeros(round(1.0 * fs), 1);
        play_and_record(hw, long_silence, spk, cfg);

        % 容器
        irReps = zeros(Lh, numErrMics, cfg.repetitions);
        allSNRs = zeros(cfg.repetitions, numErrMics);
        allCohs = zeros(cfg.repetitions, numErrMics);

        % ========== 多次重复采集 ==========
        for rep = 1:cfg.repetitions
            if rep > 1
                fprintf('  重置缓冲区...\n');
                play_and_record(hw, long_silence, spk, cfg);
                pause(0.5);
            end

            % 播放并录制
            recordedFull = play_and_record(hw, spkSig, spk, cfg);
            y_rec = recordedFull(:, errMicIdx); % [N x numErrMics]

            % >>>>>> 使用第一个通道反卷积定位主峰 <<<<<<
            inv_filter = fliplr(conj(exciteInfo.sweepSig));
            % inv_filtRer = flipud(conj(exciteInfo.sweepSig(:))); % use flipud for column vector
            Nfft = 2^nextpow2(length(recordedFull) + length(exciteInfo.sweepSig) - 1);
            
            % 对第一个麦克风做快速反卷积以定位主峰
            H_temp = fft(y_rec(:,1), Nfft) .* fft(inv_filter, Nfft);
            h_temp = ifft(H_temp, Nfft, 'symmetric');
            h_temp = h_temp(1:length(recordedFull)); % 截断至 recordedFull 长度

            % 在合理物理延迟范围内搜索主峰（10ms ～ 2.0s）
            min_delay_samples = round(0.01 * fs);   % 最小 10ms
            max_delay_samples = min(round(2.0 * fs), length(h_temp)); % 最大 2s 或信号长度
            if max_delay_samples <= min_delay_samples
                max_delay_samples = length(h_temp);
            end

            [~, peak_local_idx] = max(abs(h_temp(min_delay_samples:max_delay_samples)));
            peak_idx = min_delay_samples + peak_local_idx - 1; % 绝对索引 in recordedFull

            % 裁剪窗口：确保包含完整响应
            half_win = round(1.5 * fs);
            start_keep = max(1, peak_idx - half_win);
            end_keep = min(size(recordedFull,1), peak_idx + half_win + length(exciteInfo.sweepCore));
            y_rec_crop = y_rec(start_keep:end_keep, :);

            % 安全检查：确保裁剪后信号足够长
            if size(y_rec_crop, 1) < length(exciteInfo.sweepSig)
                warning('裁剪后信号过短 (len=%d)，可能影响IR估计', size(y_rec_crop,1));
            end

            % >>>>>> 调用 IR 估计 <<<<<<
            [irRep, irOffsets] = estimate_ir_from_sweep(...
                y_rec_crop, ...
                exciteInfo.sweepSig, ...   % 完整激励（含静音）
                exciteInfo.sweepCore, ...  % 核心扫频（预留）
                cfg);
            
            irOffsets = irOffsets + start_keep; % 转换为全局索引

            fprintf('  IR offsets(反卷积结果中的索引): [%d, %d]\n', irOffsets(1), irOffsets(2));

            irReps(:, :, rep) = irRep;

            % --- SNR与Coherence评估 ---
            snrRep = estimate_snr(y_rec, exciteInfo.sweepSig, fs);
            % ✅ 统一调用：使用 coreSweep + IR 预测
            cohRep = compute_coherence(y_rec, exciteInfo.sweepCore, irRep, irOffsets, cfg, 'Verbose', true);

            % 如果当前 rep 的平均 SNR 极低则认为该 rep 失败
            if all(snrRep < cfg.minValidSNR)
                fprintf('  ⚠️ Rep %d: SNR极低 ([%s] dB)，标记为无效\n',rep, strjoin(compose('%.1f', snrRep), ', '));
                allSNRs(rep, :) = NaN;
                allCohs(rep, :) = NaN;
            else
                allSNRs(rep, :) = snrRep;
                allCohs(rep, :) = cohRep;
            end

            fprintf('  Rep %d: SNR=[%s] dB, Coh=[%s]\n', ...
                rep, ...
                strjoin(compose('%.1f', snrRep), ', '), ...
                strjoin(compose('%.3f', cohRep), ', '));

            if rep < cfg.repetitions
                pause(cfg.repPauseSec);
            end
        end

        % ========== IR 平均（基于 SNR 加权）==========
        irAvg = robust_average_ir(irReps, allSNRs);
        irAvg = postprocess_ir(irAvg, cfg); 

        % ========== 质量综合评估 ==========
        metrics = assess_ir_quality(irAvg, allSNRs, allCohs, cfg);
        impulseResponses(:, :, spk) = irAvg;

        if metrics.usable
            fprintf('  ✅ 可用 | 中值SNR=%.1f dB, 平均Coherence=%.3f\n', ...
                metrics.medianSNR, metrics.meanCoherence);
        else
            fprintf('  ⚠️ 警告！SNR=%.1f dB, Coherence=%.3f (阈值: SNR≥%.0f, Coh≥%.2f)\n', ...
                metrics.medianSNR, metrics.meanCoherence, ...
                cfg.snrThresholdDB, cfg.coherenceThreshold);
            
            if metrics.medianSNR >= cfg.snrThresholdDB && metrics.meanCoherence < cfg.coherenceThreshold
                fprintf('  ℹ️  SNR良好但相干性低，可能是对齐问题\n');
            end
        end

        allSNRs_all(:, :, spk) = allSNRs;
        allCohs_all(:, :, spk) = allCohs;

        if spk < cfg.numSpeakers
            pause(cfg.spkSwitchPauseSec);
        end
    end

    %% ============== 保存结果 ==============
    save_secondary_results(impulseResponses, cfg, exciteInfo, allSNRs_all, allCohs_all);

catch ME
    handle_measurement_error(ME, global_hw, impulseResponses, fs);
end

fprintf('\n[complete] 测量结束。\n');
end

%% ============== 辅助函数 ==============
function cfg = validate_config_final(cfg)
% validate_config_final - 最终验证配置
    requiredFields = {'fs', 'numSpeakers', 'irMaxLen', 'micChannels', 'repetitions'};
    for f = requiredFields
        if ~isfield(cfg, f{1})
            error('缺少配置字段: %s', f{1});
        end
    end
    
    if ~isfield(cfg, 'micChannels') || ~isfield(cfg.micChannels, 'error')
        error('需指定 cfg.micChannels.error');
    end
    
    % 设置默认值
    defaults = {
        'snrThresholdDB', 12;
        'coherenceThreshold', 0.75;
        'deconvPreDelayKeep', 128;
        'repPauseSec', 1.0;
        'spkSwitchPauseSec', 2.0;
        'sweepMaxAmp', 0.9;
        'sweepMinAmp', 0.1;
        'spkAmplitude',0.8;
        'secondaryPathFile', 'secondary_path_measurement.mat'
    };
    
    for i = 1:size(defaults, 1)
        if ~isfield(cfg, defaults{i,1})
            cfg.(defaults{i,1}) = defaults{i,2};
        end
    end
end

function handle_measurement_error(ME, hw, impulseResponses, fs)
% handle_measurement_error - 处理测量错误
    fprintf('\n[FATAL ERROR] %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (%d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    
    % 尝试保存部分结果
    try
        if exist('impulseResponses', 'var') && ~isempty(impulseResponses) && ~all(impulseResponses(:) == 0)
            partialSaveFile = sprintf('secondary_path_partial_%s.mat', ...
                datetime('now', 'TimeZone', 'local'));
            save(partialSaveFile, 'impulseResponses', 'fs', '-v7.3');
            fprintf('[info] 部分结果已保存至: %s\n', partialSaveFile);
        end
    catch
    end
    
    % 释放硬件资源
    if ~isempty(hw)
        try
            hw.release();
            fprintf('[hardware-measure] 设备已释放\n');
        catch
        end
    end
    
    rethrow(ME);
end

function ir_post = postprocess_ir(ir, ~)
% postprocess_ir - 简单的后处理
% 简单后处理：窗截断
L = size(ir, 1);
win = tukeywin(L, 0.1)'; % 平滑尾部
ir_post = ir .* win';
end

function sig_norm = safe_normalize(sig, cfg)
% safe_normalize - 安全地归一化信号幅度
    maxAmp = cfg.sweepMaxAmp;
    minAmp = cfg.sweepMinAmp;
    
    sig_norm = sig(:);
    amp = max(abs(sig_norm));
    
    if amp == 0
        return;
    elseif amp > maxAmp
        sig_norm = sig_norm * (maxAmp / amp);
    elseif amp < minAmp
        sig_norm = sig_norm * (minAmp / amp);
    end
end