function [avgFreqHz, maxFreqHz, secondsAnalyzed] = analyze_wav_freq(wavPath, analyzeSeconds, method)
% per_second_avg_and_max_freq
% 逐秒计算前 N 秒音频的平均频率（谱质心）与最大频率（主峰）。
%
% 用法：
%   per_second_avg_and_max_freq('audio.wav')                  % 默认：前40秒，welch
%   per_second_avg_and_max_freq('audio.wav', 40)              % 指定时长
%   per_second_avg_and_max_freq('audio.wav', 30, 'fft')       % 指定时长与方法（fft 或 welch）
%
% 返回：
%   avgFreqHz         : 每秒平均频率（谱质心），列向量（Hz）；NaN 表示该秒能量过低或无有效谱
%   maxFreqHz         : 每秒最大频率（主峰），列向量（Hz）；NaN 表示该秒能量过低或无有效主峰
%   secondsAnalyzed   : 实际分析的秒数（可能小于请求时长）
%
% 说明：
%   - 平均频率采用谱质心：avg = sum(f .* S) / sum(S)，S 为非负谱（Welch 的 PSD 或 FFT 幅度）。
%   - 最大频率为该秒谱的最大值所在频率（排除 0 Hz）。
%   - 建议使用 'welch' 以获得更稳健结果；'fft' 更快但更易受泄漏影响。
%   - 对静音或能量很低的片段，返回 NaN；可调 minEnergyNorm 改变判断门限。

    if nargin < 2 || isempty(analyzeSeconds), analyzeSeconds = 40; end
    if nargin < 3 || isempty(method),         method = 'welch';     end

    if ~isfile(wavPath), error('文件不存在：%s', wavPath); end
    [y_full, fs] = audioread(wavPath);

    % 仅取前 analyzeSeconds 秒（若更短则截断）
    maxSamples = min(length(y_full), floor(analyzeSeconds * fs));
    y = y_full(1:maxSamples, :);
    secondsAnalyzed = floor(maxSamples / fs);

    % 多通道转单通道
    if size(y, 2) > 1
        y = mean(y, 2);
    end

    % 去直流
    y = y - mean(y);

    % 若不足1秒，直接返回空
    if secondsAnalyzed < 1
        avgFreqHz = [];
        maxFreqHz = [];
        warning('音频长度不足1秒，无法逐秒计算。');
        return;
    end

    avgFreqHz = nan(secondsAnalyzed, 1);
    maxFreqHz = nan(secondsAnalyzed, 1);
    secLen = fs;

    % 最低能量门限（基于归一化谱的最大值），避免静音秒误报
    minEnergyNorm = 0.02; % 可调整；设为 0 则不使用门限

    for s = 1:secondsAnalyzed
        idxStart = (s-1)*secLen + 1;
        idxEnd   = s*secLen;
        seg = y(idxStart:idxEnd);

        switch lower(method)
            case 'fft'
                [freqs, spec] = spectrum_fft_nonneg(seg, fs);
            case 'welch'
                [freqs, spec] = spectrum_welch_nonneg(seg, fs);
            otherwise
                error('不支持的方法：%s（可选：fft, welch）', method);
        end

        % 排除 0 Hz
        mask = freqs > 0;
        freqs = freqs(mask);
        spec  = spec(mask);

        if isempty(freqs)
            avgFreqHz(s) = NaN;
            maxFreqHz(s) = NaN;
            continue;
        end

        % 能量归一化峰值，用于最低能量判断
        peakValNorm = max(spec) / max(spec + eps); % spec 已非负，但此处仍做保护
        if peakValNorm < minEnergyNorm || sum(spec) <= eps
            avgFreqHz(s) = NaN;
            maxFreqHz(s) = NaN;
            continue;
        end

        % 最大频率（主峰）
        [~, peakIdx] = max(spec);
        maxFreqHz(s) = freqs(peakIdx);

        % 平均频率（谱质心）
        avgFreqHz(s) = sum(freqs .* spec) / sum(spec);
    end

    % 打印每秒结果
    for s = 1:secondsAnalyzed
        if isnan(avgFreqHz(s)) && isnan(maxFreqHz(s))
            fprintf('第 %2d 秒：无有效频率（能量过低或静音）\n', s);
        else
            fprintf('第 %2d 秒：平均频率 %.2f Hz，最大频率 %.2f Hz\n', s, avgFreqHz(s), maxFreqHz(s));
        end
    end
end

function [freqs, mag] = spectrum_fft_nonneg(y, fs)
% 对一秒数据做 FFT，返回非负谱（幅度）与频率轴
    N = length(y);
    if N == 0
        freqs = [];
        mag = [];
        return;
    end
    w = hann(N);
    Y = fft(y .* w);
    Nhalf = floor(N/2) + 1;
    freqs = (0:Nhalf-1)' * (fs / N);
    mag = abs(Y(1:Nhalf));      % 非负幅度谱
    % 可选：归一化到 [0,1]，但这里保留原幅度用于更稳定的质心
    % 如需归一化，请取消注释：
    % mag = mag / max(mag + eps);
end

function [freqs, Pxx] = spectrum_welch_nonneg(y, fs)
% 对一秒数据做 Welch PSD，返回非负谱（功率密度）与频率轴
    N = length(y);
    if N == 0
        freqs = [];
        Pxx = [];
        return;
    end
    winLen = max(256, round(0.5*fs)); % 0.5秒窗
    if winLen > N
        winLen = max(128, floor(N/2));
    end
    noverlap = round(0.5 * winLen);
    nfft = 2^nextpow2(winLen);

    [Pxx, freqs] = pwelch(y, hann(winLen, 'periodic'), noverlap, nfft, fs, 'onesided');
    % 返回非归一化 PSD，以获得更稳定的质心；如需归一化，可启用：
    % Pxx = Pxx / max(Pxx + eps);
end