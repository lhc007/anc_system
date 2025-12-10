function [S_hat, st_sp] = online_sec_path_estimator(st_sp, probe_sig, err_sig, cfg)
% 在线次级路径辨识（向量化 LMS 实现）
% 输入：
%   st_sp       状态结构体（首次调用时空）
%   probe_sig   [N x ns] 探测信号
%   err_sig     [N x ne] 误差信号
%   cfg         配置（含 irTruncateLen）
% 输出：
%   S_hat       [Ls x ne x ns] 更新后的次级路径估计
%   st_sp       更新后的状态

if nargin < 4 || isempty(st_sp)
    Ls = cfg.irTruncateLen;
    ne = size(err_sig,2);
    ns = size(probe_sig,2);
    st_sp.S_hat = zeros(Ls, ne, ns);
    st_sp.mu_sp = 1e-4; % 固定步长（可配置化）
end

Ls = size(st_sp.S_hat,1);
ne = size(err_sig,2);
ns = size(probe_sig,2);
N  = size(probe_sig,1);

% 向量化 LMS 更新（避免逐样本循环）
for s = 1:ns
    x = probe_sig(:,s);
    if max(abs(x)) < 1e-8, continue; end
    
    % 构造输入矩阵 X_mat [N x Ls]（因果历史）
    X_mat = zeros(N, Ls);
    for i = 1:N
        start_idx = max(1, i - Ls + 1);
        len = i - start_idx + 1;
        X_mat(i, end-len+1:end) = x(start_idx:i).';
    end
    
    for e = 1:ne
        d = err_sig(:,e);
        y_est = X_mat * st_sp.S_hat(:,e,s);
        e_lms = d - y_est;
        grad_vec = X_mat.' * e_lms; % 批量梯度
        st_sp.S_hat(:,e,s) = st_sp.S_hat(:,e,s) + st_sp.mu_sp * grad_vec;
    end
end

S_hat = st_sp.S_hat;
end