function [H, delay_samples, cfg] = load_secondary_path(varargin)
% LOAD_SECONDARY_PATH
% 加载并预处理测量得到的次级路径冲激响应
%
% 用法:
%   [H, delay_samples] = load_secondary_path;
%   [H, delay_samples] = load_secondary_path('mic', 2);       % 使用第2个麦克风
%   [H, delay_samples] = load_secondary_path('margin', 64);   % 额外延迟裕量
%
% 输出:
%   H              : [L x Ns] 复数或实数矩阵，Ns=扬声器数，L=滤波器长度
%   delay_samples  : 实际使用的对齐延迟（样本数）
%   cfg            : 原始配置结构体（含 fs 等）

p = inputParser;
addParameter(p, 'mic', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'margin', 32, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parse(p, varargin{:});

micIdx = p.Results.mic;
margin = p.Results.margin;

% --- 加载数据 ---
try
    data = load('secondary_path_measured.mat');
    if ~isfield(data, 'secondary') || ~isfield(data, 'cfg')
        error('文件 secondary_path_measured.mat 缺少必要字段');
    end
catch ME
    error('无法加载 secondary_path_measured.mat: %s', ME.message);
end

secondary = data.secondary;
cfg = data.cfg;

Ns = size(secondary.impulseResponses, 3);  % 扬声器数量
Nm = size(secondary.impulseResponses, 2);  % 麦克风数量

if micIdx > Nm
    warning('指定麦克风 %d 超出范围（共 %d 个），自动使用 mic 1', micIdx, Nm);
    micIdx = 1;
end

% --- 获取推荐延迟 ---
if isfield(secondary, 'delayEstimateSamples')
    estDelays = secondary.delayEstimateSamples(:)';
else
    % 回退：使用 summary 中的 relDelay
    estDelays = [];
    for k = 1:Ns
        estDelays(k) = secondary.summary(k).relDelay;
    end
end

% --- 统一对齐到最大延迟 + margin ---
maxDelay = max(estDelays);
totalDelay = maxDelay + margin;

% 确保不超过原始 IR 长度
origLen = size(secondary.impulseResponses, 1);
if totalDelay > origLen
    warning('请求的延迟 %d 超过原始 IR 长度 %d，截断至 %d', totalDelay, origLen, origLen);
    totalDelay = origLen;
end

% --- 提取并补零对齐 ---
H = zeros(totalDelay, Ns);
for spk = 1:Ns
    ir_raw = squeeze(secondary.impulseResponses(:, micIdx, spk));
    
    % 截取前 totalDelay 样本（若不足则补零）
    ir_padded = zeros(totalDelay, 1);
    L_use = min(length(ir_raw), totalDelay);
    ir_padded(1:L_use) = ir_raw(1:L_use);
    
    H(:, spk) = ir_padded;
end

delay_samples = totalDelay;

% --- 可选：归一化能量（按需启用）---
% energy = sqrt(sum(H.^2, 1));
% H = H ./ max(energy, eps);  % 避免除零

fprintf('[load_secondary_path] 成功加载 %d 扬声器 × %d 样本 次级路径模型\n', Ns, totalDelay);
fprintf('  使用麦克风: %d / %d\n', micIdx, Nm);
fprintf('  对齐延迟: %d 样本 (原最大估计: %d, 裕量: %d)\n', totalDelay, maxDelay, margin);
fprintf('  采样率: %d Hz → 延迟 ≈ %.1f ms\n', cfg.fs, 1000 * totalDelay / cfg.fs);

end