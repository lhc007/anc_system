function [irTrunc, energyRatio] = anc_energy_truncate(ir, keepRatio, minLen)
% ANC_ENERGY_TRUNCATE 按累积能量比例截断冲激响应
%
% 输入:
%   ir        : 列向量，冲激响应
%   keepRatio : 要保留的能量比例（如 0.99 表示保留 99% 能量）
%   minLen    : 最小输出长度（防止过短）
%
% 输出:
%   irTrunc      : 截断后的 IR
%   energyRatio  : 实际保留的能量比例

if nargin < 3, minLen = 64; end
if nargin < 2, keepRatio = 0.99; end

ir = ir(:);
L = length(ir);

if L == 0
    irTrunc = [];
    energyRatio = 0;
    return;
end

% 计算能量包络
energyEnv = abs(ir).^2;
totalEnergy = sum(energyEnv);
if totalEnergy == 0
    irTrunc = ir(1:min(minLen, L));
    energyRatio = 0;
    return;
end

% 累积能量
cumEnergy = cumsum(energyEnv);
targetEnergy = keepRatio * totalEnergy;

% 找到满足累积能量的位置
idx = find(cumEnergy >= targetEnergy, 1, 'first');

if isempty(idx)
    idx = L;
end

% 至少保留 minLen
idx = max(idx, minLen);
idx = min(idx, L);

irTrunc = ir(1:idx);
energyRatio = cumEnergy(idx) / totalEnergy;

end