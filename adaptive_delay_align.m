function best_delay = adaptive_delay_align(W, grad_buffer, cfg)
% 自适应延迟对齐
% 基于权重梯度能量动态选择最佳起始 tap
% 输入：
%   W: 当前权重 [L x nr x ns]
%   grad_buffer: 最近 K 帧的梯度累积 [L x nr x ns]
%   cfg: 配置
% 输出：
%   best_delay: 建议的 delayD（相对于滤波器末尾）

if ~isfield(cfg, 'delaySearchRange')
    cfg.delaySearchRange = [64, 256]; % 搜索范围（tap 索引）
end

L = size(W,1);
energy = sum(grad_buffer(:).^2, 3); % [L x nr]
energy_smooth = movmean(energy, [0, 32], 1); % 向后平滑

% 在指定范围内找最大能量位置
search_idx = cfg.delaySearchRange(1):min(cfg.delaySearchRange(2), L);
[~, best_tap] = max(energy_smooth(search_idx, :), [], 'all');

% 转换为 delayD：从该 tap 到末尾的样本数
best_delay = L - search_idx(best_tap);

% 限制合理范围
best_delay = max(min(best_delay, cfg.maxDelaySamples), cfg.minDelaySamples);
end