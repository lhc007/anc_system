function [y_out, st, dbg] = fxlms_recursive_process(st, refFrame, errFrame, cfg)
% 高性能帧级 FXLMS 处理（与 init 兼容）
% 主要增强：严格维护 st.xfHistory 跨帧历史，确保梯度包含跨帧贡献
% 输入:
%   st: 状态结构体（含 W{r,s,e}, S{e,s}, zi{r,s,e} 可选）
%   refFrame: [N x nr]
%   errFrame: [N x ne]
%   cfg: 配置（numRefs,numSpeakers,numErrs,timeFilterLen,useNLMS,muMax,weight_decay,maxWeightAbs,fxlmsConvMode）
% 输出:
%   y_out: [N x ns]
%   st: 更新后的状态
%   dbg: 调试信息 (W_norm, W_maxAbs, muUsed)

% 参数与校验
nr = cfg.numRefs;
ns = cfg.numSpeakers;
ne = cfg.numErrs;
N  = size(refFrame,1);
Lw = cfg.timeFilterLen;

% 输出预分配
y_out = zeros(N, ns);

% 确保 W 存在
if ~isfield(st, 'W') || isempty(st.W)
    st.W = cell(nr, ns, ne);
    for r=1:nr, for s=1:ns, for e=1:ne, st.W{r,s,e} = zeros(Lw,1); end; end; end
end

% 确认 S 存在
if ~isfield(st, 'S') || isempty(st.S)
    error('st.S 未定义，请先运行 fxlms_recursive_init');
end
useCellS = iscell(st.S);

% 选择计算模式
convMode = getfieldOrDefault(cfg, 'fxlmsConvMode', 'filter'); % 默认用 filter 模式以保持连续性

% 初始化或确认 xfHistory（保存每个 (r,s,e) 最后 Lw-1 个 Xf 样本）
if Lw > 1
    if ~isfield(st, 'xfHistory') || isempty(st.xfHistory)
        st.xfHistory = zeros(Lw-1, nr, ns, ne);
    else
        % 兼容形状
        [h1, h2, h3, h4] = size(st.xfHistory);
        if h1 ~= (Lw-1) || h2~=nr || h3~=ns || h4~=ne
            st.xfHistory = zeros(Lw-1, nr, ns, ne);
        end
    end
end

% filter 模式下确保 st.zi 存在
if strcmp(convMode, 'filter')
    if ~isfield(st, 'zi') || isempty(st.zi)
        st.zi = cell(nr, ns, ne);
        for r=1:nr, for s=1:ns, for e=1:ne
            if useCellS
                ir = st.S{e,s};
            else
                ir = squeeze(st.S(:,e,s));
            end
            st.zi{r,s,e} = zeros(max(0,length(ir)-1),1);
        end; end; end
    end
end

% 统计
weightNorm_sq_acc = zeros(ns,1);
mu_acc = 0;
W_maxAbs = 0;

% 主循环：按 (s,e) 批处理所有参考 r
for s = 1:ns
    for e_idx = 1:ne
        % 计算当前 (s,e) 的所有 nr 条 filtered-x 序列 Xf_cols [N x nr]
        Xf_cols = zeros(N, nr);
        for r = 1:nr
            if useCellS
                ir_se = st.S{e_idx, s}(:);
            else
                ir_se = squeeze(st.S(:, e_idx, s));
            end

            if strcmp(convMode, 'filter')
                zi = st.zi{r,s,e_idx};
                [xf_sig, zi_new] = filter(ir_se(:).', 1, refFrame(:, r), zi);
                Xf_cols(:, r) = xf_sig;
                st.zi{r,s,e_idx} = zi_new;
            else
                xf_full = conv(refFrame(:, r), ir_se, 'full');
                Xf_cols(:, r) = xf_full(1:min(N, length(xf_full)));
            end
        end

        % 误差信号
        d = errFrame(:, e_idx); % [N x 1]

        % 对每个参考通道 r，构造 xf_padded = [history; current] 并计算梯度（跨帧）
        for r = 1:nr
            xf_curr = Xf_cols(:, r); % [N x 1]

            if Lw > 1
                hist = squeeze(st.xfHistory(:, r, s, e_idx)); % [Lw-1 x 1]
                xf_padded = [hist(:); xf_curr(:)]; % length Lw-1+N
            else
                xf_padded = xf_curr;
            end

            % 输出贡献：conv(xf_curr, w) 只需用当前帧的 xf_curr（因果输出）
            w = st.W{r, s, e_idx};
            y_full = conv(xf_curr, w, 'full');
            y_out(:, s) = y_out(:, s) + y_full(1:N);

            % 梯度： grad[k] = sum_{n=1..N} d(n) * xf_padded((Lw-1)+n - (k-1)), k=1..Lw
            if Lw > 0
                % 使用 xcorr(xf_padded, d, Lw-1) 并按公式取 grad = r(end:-1:Lw)
                r_xf_d = xcorr(xf_padded, d, Lw-1); % lags -Lw+1..Lw-1
                % 索引说明： lag 0 -> index Lw ; lag Lw-1 -> index 2*Lw-1
                % grad for k=0..Lw-1 is r at lags (Lw-1 - k) -> indices (2Lw-1 - k)
                grad = r_xf_d((2*Lw-1) : -1 : Lw); % yields [Lw x 1] (end:-1:Lw)
                grad = grad(:);
            else
                grad = 0;
            end

            % 步长（帧级近似）
            if isfield(cfg,'useNLMS') && cfg.useNLMS
                power = sum(xf_curr.^2) + getfieldOrDefault(st,'nlmsEps',1e-8);
                mu = min(cfg.muMax, cfg.muMax / power);
            else
                mu = cfg.muMax;
            end

            % 更新权重（批量）
            w_new = w + mu * grad;

            % 衰减与裁剪
            if isfield(cfg, 'weight_decay') && cfg.weight_decay > 0
                w_new = w_new * (1 - cfg.weight_decay);
            end
            if isfield(cfg, 'maxWeightAbs') && isfinite(cfg.maxWeightAbs) && cfg.maxWeightAbs > 0
                w_new = max(min(w_new, cfg.maxWeightAbs), -cfg.maxWeightAbs);
            end

            st.W{r, s, e_idx} = w_new;

            % 统计
            weightNorm_sq_acc(s) = weightNorm_sq_acc(s) + sum(w_new.^2);
            mu_acc = mu_acc + mu;
            W_maxAbs = max(W_maxAbs, max(abs(w_new)));

            % 更新 xfHistory（保留最后 Lw-1 个 Xf 样本）
            if Lw > 1
                total_len = size(xf_padded,1);
                last_seg = xf_padded(total_len - (Lw-1) + 1 : total_len);
                st.xfHistory(:, r, s, e_idx) = last_seg;
            end
        end
    end
end

% 调试输出
dbg.W_norm = sqrt(sum(weightNorm_sq_acc) / max(1, (nr * ne * N)));
dbg.W_maxAbs = W_maxAbs;
dbg.muUsed = mu_acc / max(1, (ns * N));

end

%% 小工具
function val = getfieldOrDefault(s, field, default)
    if isstruct(s) && isfield(s, field), val = s.(field); else val = default; end
end