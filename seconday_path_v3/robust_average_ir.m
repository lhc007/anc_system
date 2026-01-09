function ir_avg = robust_average_ir(irReps, snrList)
% robust_average_ir - 基于SNR的加权平均（你已有的函数）
% irReps: [Lh x M x R]
% snrList: [R x M]
[Lh, M, R] = size(irReps);

% 每次重复的权重 = 该次所有麦克风 SNR 的平均值
weights = mean(snrList, 2); % [R x 1]
weights = weights - min(weights); % 非负化
weights = max(weights, 0);

if sum(weights) > 0
    weights = weights / sum(weights); % 归一化
else
    weights = ones(R, 1) / R; % 退化为均匀平均
end

% 初始化
ir_avg = zeros(Lh, M);
for r = 1:R
    ir_avg = ir_avg + weights(r) * irReps(:, :, r); % 直接索引，避免 squeeze
end
end