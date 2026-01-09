function [sweepSig, exciteInfo] = generate_excitation_signal(cfg)
    fs = cfg.fs;
    T_sweep = cfg.sweepDuration;            % 核心扫频时长 (秒)
    padLeading = cfg.padLeading;            % 前静音 (秒)
    padTrailing = cfg.padTrailing;          % 后静音 (秒)

    N_sweep = round(T_sweep * fs);
    N_lead = round(padLeading * fs);
    N_trail = round(padTrailing * fs);

    t = (0:N_sweep-1)' / fs;
    f_start = cfg.sweepFreqStartHz;
    f_end = cfg.sweepFreqEndHz;

    % 指数扫频（ESS）
    L = log(f_end / f_start);
    sweepCore = sin(2*pi * f_start * t .* (exp(L*t/T_sweep) - 1) / (L/T_sweep));

    % 归一化并加窗（防止边缘突变）
    win = tukeywin(N_sweep, 0.1)';
    sweepCore = sweepCore .* win';
    sweepCore = sweepCore / max(abs(sweepCore)) * cfg.amplitude;

    % 构建完整信号
    sweepSig = [zeros(N_lead, 1); sweepCore(:); zeros(N_trail, 1)];

    % === 输出信息：确保包含所有日志所需字段 ===
    exciteInfo.fs = fs;
    exciteInfo.T = T_sweep;                          % 扫频时长（秒）
    exciteInfo.coreLength = N_sweep;                 % 核心扫频样本数
    exciteInfo.preSilence = N_lead;                  % ✅ 新增：前导静音样本数
    exciteInfo.postSilence = N_trail;                % ✅ 新增：后导静音样本数
    exciteInfo.sweepStartIdx = N_lead + 1;           % 扫频起始索引（1-based）
    exciteInfo.totalLength = length(sweepSig);       % 总长度（样本）
    exciteInfo.sweepCore = sweepCore(:);             % 核心扫频信号
    exciteInfo.sweepSig = sweepSig;                  % 完整激励信号（含静音）
end