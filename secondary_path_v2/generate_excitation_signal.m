function [sweepSig, info] = generate_excitation_signal(cfg)
    fs = cfg.fs;
    T_sweep = cfg.sweepDuration;            % 核心扫频时长 (秒)
    padLeading = cfg.padLeading;            % 前静音 (秒)
    padTrailing = cfg.padTrailing;          % 后静音 (秒)

    N_sweep = round(T_sweep * fs);
    N_lead = round(padLeading * fs);
    N_trail = round(padTrailing * fs);

    t = (0:N_sweep-1)' / fs;
    f_start = cfg.f_min;
    f_end = cfg.f_max;

    % 指数扫频（ESS）
    L = log(f_end / f_start);
    sweepCore = sin(2*pi * f_start * t .* (exp(L*t/T_sweep) - 1) / (L/T_sweep));

    % 归一化并加窗（防止边缘突变）
    win = tukeywin(N_sweep, 0.1)';
    sweepCore = sweepCore .* win';
    sweepCore = sweepCore / max(abs(sweepCore)) * cfg.amplitude;

    % 构建完整信号
    sweepSig = [zeros(N_lead, 1); sweepCore(:); zeros(N_trail, 1)];

    % 输出信息
    info.fs = fs;
    info.coreLength = N_sweep;
    info.sweepStartIdx = N_lead + 1;
    info.sweepSignal = sweepCore(:);  % 仅核心部分
    info.totalLength = length(sweepSig);
end