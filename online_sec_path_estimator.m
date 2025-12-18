function [S_hat, st_sp] = online_sec_path_estimator(st_sp, probe_sig, err_sig, cfg)
% 在线次级路径估计（帧级向量化 LMS）
% 使用互相关批量计算梯度，兼顾效率与稳定性
% 输入:
%   st_sp: 状态，可为空（首次调用）
%   probe_sig: [N x ns]
%   err_sig:   [N x ne]
%   cfg: 包含 irTruncateLen, mu_sp (可选)
% 输出:
%   S_hat: [Ls x ne x ns]
%   st_sp: 更新后的状态

if nargin < 4, error('online_sec_path_estimator: 需要 st_sp, probe_sig, err_sig, cfg'); end

[N, ns] = size(probe_sig);
[~, ne] = size(err_sig);

% 初始化状态
if isempty(st_sp)
    Ls = cfg.irTruncateLen;
    st_sp.S_hat = zeros(Ls, ne, ns);
    st_sp.mu_sp = getfieldOrDefault(cfg, 'mu_sp', 1e-4);
end
Ls = size(st_sp.S_hat,1);
mu_sp = getfieldOrDefault(st_sp, 'mu_sp', getfieldOrDefault(cfg, 'mu_sp', 1e-4));

% 对每个扬声器并行处理（向量化）
for s = 1:ns
    x = probe_sig(:, s);
    if max(abs(x)) < 1e-12
        continue;
    end

    % 预计算能量用于步长归一化
    energy = sum(x.^2) + 1e-12;

    for e = 1:ne
        d = err_sig(:, e);

        % 使用互相关计算批量梯度： grad_k = sum_n d(n) * x(n-k+1)
        r_x_d = xcorr(x, d, Ls-1); % length 2*Ls-1
        grad = r_x_d(Ls : Ls + Ls - 1); % [Ls x 1], lags 0..Ls-1

        % 自适应步长（依据 probe 能量）
        mu_eff = mu_sp / energy;

        % 更新
        st_sp.S_hat(:, e, s) = st_sp.S_hat(:, e, s) + mu_eff * grad;
    end
end

S_hat = st_sp.S_hat;

end

%% 小工具
function val = getfieldOrDefault(s, field, default)
    if isstruct(s) && isfield(s, field), val = s.(field); else val = default; end
end