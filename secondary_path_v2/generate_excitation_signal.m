function [sweepCore_scaled, sweepSig, info, weightedSweepCore, invFilter] = generate_excitation_signal(cfg)
%GENERATE_EXCITATION_SIGNAL 生成ESS激励信号（带同步Click）
%
% 输出：
%   sweepCore_scaled    : 缩放后的核心扫频（用于播放）
%   sweepSig            : 完整激励信号（含前后静音 + Click）
%   info                : 信号参数结构体
%   weightedSweepCore   : 加权扫频核心（用于相干性重建）
%   invFilter           : 反卷积滤波器（必须基于 scaled 信号！）

    % --- 参数解析 ---
    fs = cfg.fs;
    T = cfg.sweepDuration;          % 核心扫频时长 (秒)
    f1 = cfg.sweepFreqStartHz;      % 起始频率
    f2 = cfg.sweepFreqEndHz;        % 结束频率
    padLeading = cfg.padLeading;    % 前导静音 (秒)
    padTrailing = cfg.padTrailing;  % 尾随静音 (秒)

    sweepCoreLen = round(T * fs);
    padLeadSamples = round(padLeading * fs);
    padTrailSamples = round(padTrailing * fs);

    % --- 生成ESS核心（仅信号，不包含 invFilter）---
    [sweepCore] = generate_ess_core_signal(f1, f2, sweepCoreLen, fs);

    % --- 幅度缩放 ---
    maxAmp = cfg.amplitude;
    scale = maxAmp / max(abs(sweepCore));
    sweepCore_scaled = scale * sweepCore;

    % --- 计算加权函数（与 sweepCore 同长度）---
    t = (0:sweepCoreLen-1)' / fs;
    weight = exp(-t * log(f2/f1) / T);

    % weightedSweepCore 用 scaled 版本（用于相干性）
    weightedSweepCore = weight .* sweepCore_scaled;

    % === 关键修复：invFilter 必须基于实际播放的 scaled 信号！===
    invFilter = flipud(weightedSweepCore);  % ← 正确做法！

    % --- 构建完整激励信号 ---
    totalLen = padLeadSamples + sweepCoreLen + padTrailSamples;
    sweepSig = zeros(totalLen, 1);

    % === 插入同步 Click ===
    click_offset_sec = cfg.click_offset_sec;
    click_pos = padLeadSamples - round(click_offset_sec * fs); % 在前静音末尾
    if click_pos > 0 && click_pos < padLeadSamples
        sweepSig(click_pos) = 0.8; % 幅度 0.8，避免削波
        fprintf('[signal] 已插入同步 Click @ %.3f 秒\n', click_pos/fs);
    end

    % === 放置核心扫频 ===
    sweepSig(padLeadSamples + (1:sweepCoreLen)) = sweepCore_scaled;

    % --- 返回信息 ---
    info.coreLength = sweepCoreLen;
    info.totalLength = totalLen;
    info.padLeadingSamples = padLeadSamples;
    info.padTrailingSamples = padTrailSamples;
    info.fs = fs;
    info.T = T;
    info.clickPos = click_pos; % 供 align_sweep_start 使用
    info.sweepStartIdx = padLeadSamples + 1; % 扫频理论起始索引
    info.sweepStartSec = padLeading; % 扫频理论起始时间（秒）
end

% --- 辅助函数：仅生成 ESS 信号（无 invFilter）---
function sweepCore = generate_ess_core_signal(f1, f2, N, fs)
    if N <= 0
        error('无效的 sweepCore 长度');
    end
    t = (0:N-1)' / fs;
    phase = 2 * pi * f1 * N / (fs * log(f2/f1)) * (exp(t * log(f2/f1) * fs / N) - 1);
    % 更稳健的相位计算（避免 T=0 问题）
    sweepCore = sin(phase);
end