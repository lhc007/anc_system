function [irAligned, delayApplied, peakIdx] = anc_path_alignment(ir, threshDb, offset, targetIdx, maxLength)
% ANC_PATH_ALIGNMENT 对冲激响应进行峰值检测与对齐
%
% 输入:
%   ir          : 列向量，原始冲激响应
%   threshDb    : 峰值检测阈值（dB，相对于最大值）
%   offset      : 在峰值前额外保留的样本数（用于因果性）
%   targetIdx   : 目标对齐位置（如 256），若 <=0 则自动计算
%   maxLength   : 输出最大长度
%
% 输出:
%   irAligned     : 对齐后的 IR（长度 = maxLength）
%   delayApplied  : 实际插入的零延迟（样本数）
%   peakIdx       : 原始 IR 中主峰位置（1-based）

if nargin < 5, maxLength = length(ir); end
if nargin < 4, targetIdx = -1; end
if nargin < 3, offset = 64; end
if nargin < 2, threshDb = -20; end

ir = ir(:);
L = length(ir);

% --- 峰值检测 ---
env = abs(ir);
maxVal = max(env);
thresh = maxVal * 10^(threshDb/20);
valid = env >= thresh;

if ~any(valid)
    % 无有效峰值，使用最大值位置
    [~, peakIdx] = max(env);
else
    % 找第一个超过阈值的点（避免尾部噪声）
    peakIdx = find(valid, 1, 'first');
end

% --- 计算目标起始位置 ---
if targetIdx > 0
    desiredStart = targetIdx;
else
    % 自动：将峰值放在 offset 之后
    desiredStart = peakIdx + offset;
end

% 实际需要在前面补零的数量
delayApplied = max(0, desiredStart - peakIdx);

% 构建对齐信号
irPadded = [zeros(delayApplied, 1); ir];

% 截断或补零到 maxLength
if length(irPadded) > maxLength
    irAligned = irPadded(1:maxLength);
else
    irAligned = [irPadded; zeros(maxLength - length(irPadded), 1)];
end

end