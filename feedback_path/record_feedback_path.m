function record_feedback_path()
% 录制多通道反馈路径（扬声器 → 参考麦克风）
% 使用 hardware_init_measure（流式帧接口）
% 输出: feedback_path.mat (含 impulseResponses_fb 和 delayEstimate_fb)
%
% ⚠️ 必须在 ANC 关闭状态下运行！

    cfg = anc_config(); % 加载配置
    
    %% === 参数提取 ===
    fs           = cfg.fs;
    numSpk       = cfg.numSpeakers;               % 应为 4
    refCh        = cfg.micChannels.reference;     % e.g., [1 2 3 4]
    numRef       = numel(refCh);
    
    sweepDur     = cfg.fbSweepDur;                % 扫频时长 (秒)
    fStart       = cfg.fbSweepFstart;             % 起始频率 (Hz)
    fEnd         = cfg.fbSweepFend;               % 截止频率 (Hz)
    nSweeps      = cfg.fbSweepNsweeps;            % 扫频次数
    useFreqWgt   = cfg.fbUseFreqWeight;           % 是否低频加权
    lowBoostHz   = cfg.fbLowBoostHz;
    lowBoostPwr  = cfg.fbLowBoostPower;
    
    alignTarget  = cfg.fbAlignTargetIdx;
    alignOffset  = cfg.fbAlignOffset;
    alignThresh  = cfg.fbAlignThreshDb;
    minCorr      = cfg.fbCorrMin;
    energyRatio  = cfg.fbEnergyCutRatio;
    irMaxLen     = cfg.irMaxLen;
    irTruncLen   = cfg.irTruncateLen;

    %% === 安全生成扫频信号（频域加权，避免 Inf）===
    t = (0:1/fs:sweepDur - 1/fs).';
    sweepSig = chirp(t, fStart, sweepDur, fEnd, 'logarithmic');
    
    if useFreqWgt
        N = length(sweepSig);
        f_fft = (0:N-1)' * fs / N;
        f_fft = f_fft(1:floor(N/2)+1);  % 单边频谱
        
        weight = ones(size(f_fft));
        idx = (f_fft > 0) & (f_fft <= lowBoostHz);  % 跳过 f=0
        weight(idx) = (f_fft(idx) / lowBoostHz).^(-lowBoostPwr);
        
        % 对称扩展至全频（用于时域滤波）
        if mod(N, 2) == 0
            weight_full = [weight; flipud(weight(2:end-1))];
        else
            weight_full = [weight; flipud(weight(2:end))];
        end
        
        sweep_weighted = ifft(fft(sweepSig) .* weight_full);
        sweepSig = real(sweep_weighted);
    end
    
    sweepSig = sweepSig / max(abs(sweepSig) + 1e-12); % 防除零
    L_sweep = length(sweepSig);
    
    % ✅ 强制检查 finite
    if ~all(isfinite(sweepSig))
        error('[feedback] sweepSig contains NaN or Inf!');
    end

    %% === 初始化硬件（流式帧接口）===
    fprintf('🔧 初始化硬件...\n');
    hw = hardware_init_measure(cfg);
    
    % 清空设备缓冲区
    for i = 1:5
        hw.reader();
    end
    pause(0.1);

    %% === 录制参数 ===
    extraTime = 0.2; % 额外录制时间（捕获 IR 尾部）
    totalRecSamples = ceil((sweepDur + extraTime) * fs);
    frameSize = cfg.timeFrameSamples;
    allIR = zeros(irMaxLen, numRef, numSpk, nSweeps);

    fprintf('\n🔁 开始录制反馈路径...\n');
    fprintf(' 扬声器: %d, 参考麦: %d, 扫频长度: %.1f s\n', numSpk, numRef, sweepDur);

    %% === 主循环：扬声器 × Sweep ===
    for s = 1:numSpk
        fprintf('\n▶ 扬声器 %d/%d\n', s, numSpk);
        
        for k = 1:nSweeps
            % 构建 4 通道播放信号（仅第 s 通道有信号）
            playSignal = zeros(L_sweep, numSpk);
            playSignal(:, s) = sweepSig;
            
            recBuffer = [];
            idx = 1; % 播放指针
            
            while size(recBuffer, 1) < totalRecSamples
                % --- 准备当前播放帧 ---
                if idx <= L_sweep
                    frameOut = playSignal(idx:min(idx + frameSize - 1, L_sweep), :);
                    idx = idx + size(frameOut, 1);
                else
                    needed = totalRecSamples - size(recBuffer, 1);
                    frameOut = zeros(min(frameSize, needed), numSpk);
                end
                
                % 补零到 frameSize（确保 writer 接收固定长度）
                if size(frameOut, 1) < frameSize
                    frameOut(end+1:frameSize, :) = 0;
                end
                
                % --- 同步播放 + 录音 ---
                hw.writer(frameOut);
                recFrame = hw.reader(); % [frameSize × micNumChannels]
                
                recBuffer = [recBuffer; recFrame];
            end
            
            % 提取参考麦克风通道
            recRef = recBuffer(1:totalRecSamples, refCh);
            
            % 解卷积（使用维纳滤波）
            IR_raw = deconvwnr(recRef, sweepSig, 1e-3);
            
            % 补零/截断到 irMaxLen
            if size(IR_raw, 1) > irMaxLen
                IR_raw = IR_raw(1:irMaxLen, :);
            else
                IR_raw(end+1:irMaxLen, :) = 0;
            end
            
            allIR(:, :, s, k) = IR_raw;
            fprintf('  Sweep %d/%d 完成.\n', k, nSweeps);
            
            pause(0.1); % 避免设备过载
        end
    end

    %% === 释放硬件 ===
    hw.release();
    fprintf('✅ 硬件已释放。\n');

    %% === 后处理：对齐 + 平均 ===
    fprintf('\n🛠️  后处理反馈路径...\n');
    alignedIR = zeros(irMaxLen, numRef, numSpk);
    delayEst = zeros(numRef, numSpk);

    for s = 1:numSpk
        for r = 1:numRef
            irSet = squeeze(allIR(:, r, s, :));
            irRef = mean(irSet, 2);
            
            alignedSet = zeros(size(irSet));
            corrs = zeros(1, nSweeps);
            
            for k = 1:nSweeps
                ir_k = irSet(:, k);
                xc = xcorr(irRef, ir_k, 100, 'coeff');
                [~, lagIdx] = max(abs(xc));
                lag = lagIdx - 101;
                corrs(k) = max(abs(xc));
                
                if lag > 0
                    alignedSet(:, k) = [zeros(lag,1); ir_k(1:end-lag)];
                else
                    alignedSet(:, k) = [ir_k(1-lag:end); zeros(-lag,1)];
                end
            end
            
            validMask = corrs >= minCorr;
            if sum(validMask) == 0
                warning('Spk%d→Ref%d: 无有效 sweep，使用全部数据', s, r);
                validMask = true(1, nSweeps);
            end
            
            irAvg = mean(alignedSet(:, validMask), 2);
            [irFinal, ~, peakIdx] = anc_path_alignment(...
                irAvg, alignThresh, alignOffset, alignTarget, irMaxLen);
            
            alignedIR(:, r, s) = irFinal;
            delayEst(r, s) = peakIdx;
        end
    end

       %% === 截断（按能量或固定长度）===
    if irTruncLen > 0 && irTruncLen < irMaxLen
        alignedIR = alignedIR(1:irTruncLen, :, :);
    else
        for s = 1:numSpk
            for r = 1:numRef
                ir_orig = alignedIR(:, r, s);
                [ir_trunc, ~] = anc_energy_truncate(ir_orig, energyRatio, 32);
                Ltrunc = length(ir_trunc);
                alignedIR(1:Ltrunc, r, s) = ir_trunc;
                alignedIR(Ltrunc+1:end, r, s) = 0;
            end
        end
    end

    %% === 构造 F 并保存 ===
    L_f = round(cfg.fs * cfg.feedbackIRLenSec);
    currentLen = size(alignedIR, 1);
    
    if currentLen >= L_f
        F = alignedIR(1:L_f, :, :);
    else
        F = zeros(L_f, numRef, numSpk);
        F(1:currentLen, :, :) = alignedIR;
    end

    % 保存 F（这才是 ANC 系统需要的反馈路径）
    save(cfg.feedbackPathFile, 'F');
    
    fprintf('\n✅ 反馈路径录制完成！\n');
    fprintf(' 文件: %s\n', cfg.feedbackPathFile);
    fprintf(' 尺寸: [%d x %d x %d]\n', size(F));
    fprintf(' 估计延迟: %d 样本 (%.2f ms)\n', ...
        median(delayEst(:)), ...
        1000 * median(delayEst(:)) / fs);