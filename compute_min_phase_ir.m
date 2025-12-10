function ir_min = compute_min_phase_ir(ir, outLen, normalizeToPeak)
% 计算给定冲激响应的最小相位等效版本（优化工具）
% 输入：
%   ir                原始冲激响应 (列向量)
%   outLen            输出长度（可用于压缩或截断）
%   normalizeToPeak   是否将最小相位 IR 的峰值匹配原始峰值（true/false）
%
% 原理：
%   - 对原始 IR 做 FFT 得到频谱 H(ω)
%   - 取 log|H(ω)| 的 IFFT 得到倒谱
%   - 保留 0 点 + 正频率倒谱系数的 2 倍（最小相位倒谱）
%   - 再 FFT → exp(·) → IFFT 得到最小相位时域响应
%
% 注意：
%   - 最小相位 IR 并非真实物理次级路径，不能直接替换 Filtered-X 模型
%   - 可用于初始化控制器权重或压缩路径长度
%
% 作者：重构版

ir = ir(:);
N  = length(ir);
if nargin < 2 || isempty(outLen)
    outLen = N;
end
if nargin < 3
    normalizeToPeak = true;
end

Nfft = 2^nextpow2(N*2);
H    = fft(ir, Nfft);
mag  = abs(H)+1e-12;
logMag = log(mag);
cep  = real(ifft(logMag)); % 实倒谱

% 构造最小相位倒谱
cep_min = zeros(size(cep));
cep_min(1) = cep(1);
half = floor(Nfft/2);
cep_min(2:half) = 2*cep(2:half); % 正频率部分乘2

% 频域重建
H_min = exp(fft(cep_min));
h_min = real(ifft(H_min));

ir_min = h_min(1:outLen);

if normalizeToPeak
    p_orig = max(abs(ir))+1e-12;
    p_min  = max(abs(ir_min))+1e-12;
    ir_min = ir_min * (p_orig / p_min);
end
end