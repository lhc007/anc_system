function [y_out, st, dbg] = fxlms_recursive_process(st, refFrame, errFrame, cfg)
% 高性能 FXLMS 递推处理（数组化 W + 窗口向量化 + 可选 FFT 卷积）
%
% 输入:
%   st        : 状态结构体（含 W, S, xfHistory/zi_fxf 等）
%   refFrame  : [N × nr] 参考信号
%   errFrame  : [N × ne] 误差信号
%   cfg       : 配置结构体
%
% 输出:
%   y_out     : [N × ns] 控制输出
%   st        : 更新后的状态
%   dbg       : 调试信息

% ---------------- 参数与维度 ----------------
nr = cfg.numRefs;
ns = cfg.numSpeakers;
ne = cfg.numErrs;
N  = size(refFrame, 1);
Lw = cfg.timeFilterLen;

% 安全检查
assert(isequal(size(errFrame,1), N), 'errFrame 必须与 refFrame 帧长一致');
assert(isequal(size(refFrame,2), nr), 'refFrame 列数必须为 nr');
assert(isequal(size(errFrame,2), ne), 'errFrame 列数必须为 ne');

% 输出缓冲
y_out = zeros(N, ns);

% ---------------- Filtered-X 计算 ----------------
% 支持两种模式：'conv'（默认）与 'filter'
convMode = getfieldOrDefault(cfg, 'fxlmsConvMode', 'conv');

% 预分配 Filtered-X: [N × nr × ns × ne]
Xf = zeros(N, nr, ns, ne);

if strcmp(convMode, 'filter')
    % 使用 filter + zi 的有状态因果滤波
    % 需要 st.zi_fxf{e,s,r} 与 st.S_cell{e,s}
    for s = 1:ns
        for e_idx = 1:ne
            ir = st.S_cell{e_idx, s};   % [Ls_e,s × 1]
            for r = 1:nr
                [xf_sig, zi_new] = filter(ir, 1, refFrame(:, r), st.zi_fxf{e_idx, s, r});
                Xf(:, r, s, e_idx) = xf_sig;
                st.zi_fxf{e_idx, s, r} = zi_new;
            end
        end
    end
else
    % 使用时域卷积（适合中短 IR）
    for s = 1:ns
        for e_idx = 1:ne
            ir = st.S_cell{e_idx, s};   % [Ls_e,s × 1]
            for r = 1:nr
                xf_full = conv(refFrame(:, r), ir, 'full');
                Xf(:, r, s, e_idx) = xf_full(1:N);
            end
        end
    end
end

% ---------------- 窗口缓冲（长度 Lw）----------------
% winBuf: [Lw × nr × ns × ne]，首行是最新样本
if ~isfield(st, 'xfHistory') || isempty(st.xfHistory)
    st.xfHistory = zeros(max(Lw-1,0), nr, ns, ne);
end
winBuf = zeros(Lw, nr, ns, ne);
if Lw > 1
    winBuf(2:end, :, :, :) = st.xfHistory; % 历史填到底部
end

% ---------------- 统计与步长 ----------------
dbg = struct();
dbg.W_norm   = 0;
dbg.W_maxAbs = 0;
dbg.muUsed   = 0;

% 将 W 改为数组: [Lw × nr × ns × ne]
W = st.W;  % 数组（已在 init 中设置）

% ---------------- 主迭代：逐样本递推 ----------------
for n = 1:N
    % 推入新样本到窗口顶部
    % winBuf = [new; old(1:end-1)]
    winBuf = cat(1, reshape(Xf(n, :, :, :), [1, nr, ns, ne]), winBuf(1:Lw-1, :, :, :));

    % 计算输出 y_out(n,s)
    % 对每个 s，先在 Lw 维做加权内积，再对 r,e 累加
    % contrib_sre = sum_Lw(W .* winBuf) -> [nr × ne]
    for s = 1:ns
        % 取该扬声器的窗口与权重切片
        Hs = squeeze(winBuf(:, :, s, :)); % [Lw × nr × ne]
        Ws = squeeze(W(:, :, s, :));      % [Lw × nr × ne]
        % Lw 维相乘求和 -> [nr × ne]
        contrib = squeeze(sum(Ws .* Hs, 1)); % [nr × ne]
        % 与误差 e(n,:) 相乘并对 r,e 累加得到标量输出
        y_out(n, s) = sum(contrib * errFrame(n, :).');
    end

    % NLMS 步长（对每个 r,s,e 单独）
    if cfg.useNLMS
        % winBuf 能量（Lw 维求和）：[nr × ns × ne]
        power = squeeze(sum(winBuf.^2, 1)); % [nr × ns × ne]
        muMat = cfg.muMax ./ (power + 1e-12); % [nr × ns × ne]
    else
        muMat = ones(nr, ns, ne) * cfg.muMax;
    end

    % 权重更新（向量化）
    % dW = mu * (winBuf .* err_n)；将 err_n 扩展到最后一维 ne
    err_n = reshape(errFrame(n, :), [1, 1, 1, ne]); % [1 × 1 × 1 × ne]
    dW = bsxfun(@times, winBuf, err_n);             % [Lw × nr × ns × ne]
    % 再乘以 muMat（按 r,s,e）
    dW = bsxfun(@times, dW, reshape(muMat, [1, nr, ns, ne]));
    W = W + dW;

    % 衰减与裁剪（可选）
    if isfield(cfg, 'weight_decay') && cfg.weight_decay > 0
        W = (1 - cfg.weight_decay) * W;
    end
    if isfield(cfg, 'maxWeightAbs') && isfinite(cfg.maxWeightAbs) && cfg.maxWeightAbs > 0
        W = max(min(W, cfg.maxWeightAbs), -cfg.maxWeightAbs);
    end

    % 步长统计
    if cfg.useNLMS
        avg_mu = mean(muMat(:));
        dbg.muUsed = dbg.muUsed + avg_mu;
    else
        dbg.muUsed = dbg.muUsed + cfg.muMax;
    end
end

% 调试信息（范数与最大绝对值）
dbg.W_norm   = norm(W(:));
dbg.W_maxAbs = max(abs(W(:)));
dbg.muUsed   = dbg.muUsed / N;

% 帧尾更新历史
if Lw > 1
    st.xfHistory = winBuf(2:end, :, :, :);
end

% 写回权重
st.W = W;

end

% 安全获取字段
function val = getfieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end