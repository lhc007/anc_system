function [ir, snr_db, peakIdx, meta] = deconvolve_and_estimate_snr(rec_sig, sweep, fs, varargin)
% DECONVOLVE_AND_ESTIMATE_SNR
%   从录制信号中提取冲激响应并估算声学 SNR。
%
% 输入:
%   rec_sig      : 录制信号 (列向量)
%   sweep        : 原始扫频信号 (列向量)
%   fs           : 采样率 (Hz)
%   varargin     : 可选参数（键值对）
%       'sigWinRadius'    : 信号窗半宽（样本），默认 200
%       'noiseWinStart'   : 噪声窗起始偏移（主峰后多少样本），默认 2000
%       'noiseWinLength'  : 噪声窗长度（样本），默认 2000
%       'irMaxLenExtra'   : IR 最大长度 = length(sweep) + extra，默认 10000
%       'regEps'          : 反卷积正则化项，默认 1e-12
%
% 输出:
%   ir           : 冲激响应（截断后）
%   snr_db       : 估算的 SNR (dB)
%   peakIdx      : 主峰位置（在 ir 中的索引）
%   meta         : 结构体，包含对齐延迟、窗口位置等调试信息

p = inputParser;
addRequired(p, 'rec_sig', @(x) isvector(x));
addRequired(p, 'sweep', @(x) isvector(x));
addRequired(p, 'fs', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sigWinRadius', 200, @isscalar);
addParameter(p, 'noiseWinStart', 2000, @isscalar);
addParameter(p, 'noiseWinLength', 2000, @isscalar);
addParameter(p, 'irMaxLenExtra', 10000, @isscalar);
addParameter(p, 'regEps', 1e-12, @isscalar);
parse(p, rec_sig, sweep, fs, varargin{:});

rec_sig = rec_sig(:);
sweep = sweep(:);

L = length(sweep);
irMaxLen = L + p.Results.irMaxLenExtra;

%% === 1. 对齐：找 sweep 在录制信号中的起始位置 ===
fprintf('  - 对齐 sweep...\n');
[xc, lags] = xcorr(rec_sig, sweep);
[~, imax] = max(abs(xc));
delay_samples = lags(imax);  % delay = rec_start - sweep_start
start_idx = delay_samples + 1;
end_idx = start_idx + L - 1;

if start_idx < 1 || end_idx > length(rec_sig)
    error('对齐失败：响应超出录制范围。start=%d, end=%d, total=%d', ...
          start_idx, end_idx, length(rec_sig));
end

rec_seg = rec_sig(start_idx:end_idx);

%% === 2. 反卷积 ===
fprintf('  - 执行反卷积...\n');
Nfft = 2^nextpow2(2*L);
S = fft(sweep, Nfft);
InvFilter = conj(S) ./ (abs(S).^2 + p.Results.regEps);
Rec = fft(rec_seg, Nfft);
H = Rec .* InvFilter;
ir_full = real(ifft(H));

% 截断到合理长度
ir = ir_full(1:min(length(ir_full), irMaxLen));

%% === 3. 找主峰 ===
[~, peakIdx] = max(abs(ir));
if peakIdx < 100
    warning('主峰过早出现（<100 样点），可能包含电子串扰');
end

%% === 4. 计算 SNR ===
% 信号窗：主峰附近
sig_start = max(1, peakIdx - p.Results.sigWinRadius);
sig_end   = min(length(ir), peakIdx + p.Results.sigWinRadius);
signalPower = mean(ir(sig_start:sig_end).^2);

% 噪声窗：主峰之后的平稳段
noise_start_raw = peakIdx + p.Results.noiseWinStart;
noise_start = min(length(ir), noise_start_raw);
noise_end   = min(length(ir), noise_start + p.Results.noiseWinLength - 1);

% 如果主峰太靠后，fallback 到尾部
if noise_start >= length(ir)
    noise_start = max(1, length(ir) - p.Results.noiseWinLength);
    noise_end = length(ir);
end

noisePower = mean(ir(noise_start:noise_end).^2 + 1e-18);
snr_db = 10 * log10(signalPower / noisePower);

%% === 5. 返回元数据 ===
meta.delay_samples = delay_samples;        % sweep 起始延迟
meta.start_idx = start_idx;
meta.ir_length = length(ir);
meta.sig_win = [sig_start, sig_end];
meta.noise_win = [noise_start, noise_end];
meta.signal_power = signalPower;
meta.noise_power = noisePower;