function [y_out, st, dbg] = fxlms_recursive_process(st, refFrame, errFrame, cfg)
% 高性能 FXLMS 递推处理（严格匹配接口 + 修正维度/统计/边界问题）
%
% 输入:
%   st        : 状态结构体（含 W, S, zi 或 xfHistory）
%   refFrame  : [N × nr] 参考信号
%   errFrame  : [N × ne] 误差信号
%   cfg       : 配置结构体
%
% 输出:
%   y_out     : [N × ns] 控制输出（扬声器驱动）
%   st        : 更新后的状态
%   dbg       : 调试信息（weightNorm, muUsed）

% === 参数与维度 ===
nr = cfg.numRefs;
ns = cfg.numSpeakers;
ne = cfg.numErrs;
N  = size(refFrame, 1);      % 帧长度
Lw = cfg.timeFilterLen;      % 控制滤波器长度

% === 输出初始化 ===
y_out = zeros(N, ns);  % [N × ns]

% === Filtered-X 计算模式选项 ===
convMode = getfieldOrDefault(cfg, 'fxlmsConvMode', 'conv');  % 'conv' 或 'filter'

% === 预分配 Filtered-X 缓冲 ===
Xf = zeros(N, nr, ns, ne);

% === Step 1: 计算 Filtered-X ===
if strcmp(convMode, 'filter')
    % 模式 A: filter + zi（因果，有状态，适合长 IR）
    if ~isfield(st, 'zi_fxf') || isempty(st.zi_fxf)
        st.zi_fxf = cell(ne, ns, nr);
        for e_idx = 1:ne
            for s = 1:ns
                ir_se = st.S{e_idx, s};
                for r = 1:nr
                    st.zi_fxf{e_idx, s, r} = zeros(max(0, length(ir_se)-1), 1);
                end
            end
        end
    end

    for s = 1:ns
        for e_idx = 1:ne
            ir_se = st.S{e_idx, s};
            for r = 1:nr
                [xf_sig, zi_new] = filter(ir_se, 1, refFrame(:, r), st.zi_fxf{e_idx, s, r});
                Xf(:, r, s, e_idx) = xf_sig;
                st.zi_fxf{e_idx, s, r} = zi_new;
            end
        end
    end
else
    % 模式 B: conv（快，适合短 IR；取前 N 保持因果前缀）
    for s = 1:ns
        for e_idx = 1:ne
            ir_se = st.S{e_idx, s};
            for r = 1:nr
                xf_full = conv(refFrame(:, r), ir_se, 'full');
                Xf(1:N, r, s, e_idx) = xf_full(1:N);
            end
        end
    end
end

% === Step 2: 初始化历史缓冲（用于 Lw > 1 的 FIR 权重）===
if ~isfield(st, 'xfHistory') || isempty(st.xfHistory)
    st.xfHistory = zeros(max(Lw - 1, 0), nr, ns, ne);
end

% === 调试统计初始化 ===
total_weight_norm_sq = zeros(ns, 1);  % 每扬声器累计权重范数平方
total_mu = 0;

% === Step 3: 逐样本 FXLMS 更新 ===
for n = 1:N
    e_n = errFrame(n, :)';  % [ne × 1]
    
    for s = 1:ns
        y_val = 0;
        local_power = 0;
        
        for r = 1:nr
            for e_idx = 1:ne
                % --- 构造 Filtered-X 向量 (长度 = Lw) ---
                if Lw == 1
                    xf_vec = Xf(n, r, s, e_idx);
                else
                    if n == 1
                        xf_vec = [Xf(n, r, s, e_idx); st.xfHistory(:, r, s, e_idx)];
                    else
                        startIdx = max(1, n - Lw + 1);
                        recent_xf = Xf(startIdx:n, r, s, e_idx);
                        if length(recent_xf) < Lw
                            pad_len = Lw - length(recent_xf);
                            xf_vec = [recent_xf; st.xfHistory(1:pad_len, r, s, e_idx)];
                        else
                            xf_vec = recent_xf;
                        end
                    end
                end
                
                % --- 权重输出 ---
                w = st.W{r, s, e_idx};  % [Lw × 1]
                y_val = y_val + w.' * xf_vec;
                
                % --- 步长计算（NLMS 或固定）---
                if cfg.useNLMS
                    power = sum(xf_vec.^2) + 1e-12;
                    mu = cfg.muMax / power;
                else
                    mu = cfg.muMax;
                end
                
                % --- 梯度与权重更新 ---
                grad = e_n(e_idx) * xf_vec;
                w_new = w + mu * grad;
                
                % --- 稳定化：衰减 + 可选裁剪 ---
                if isfield(cfg, 'weight_decay') && cfg.weight_decay > 0
                    w_new = w_new * (1 - cfg.weight_decay);
                end
                if isfield(cfg, 'maxWeightAbs') && cfg.maxWeightAbs > 0
                    w_new = min(w_new,  cfg.maxWeightAbs);
                    w_new = max(w_new, -cfg.maxWeightAbs);
                end
                st.W{r, s, e_idx} = w_new;
                
                % --- 调试统计：累加权重范数平方 ---
                total_weight_norm_sq(s) = total_weight_norm_sq(s) + norm(w_new)^2;
                
                % 累加功率（用于平均 mu）
                local_power = local_power + sum(xf_vec.^2);
            end
        end
        
        % 赋值输出（[N × ns]）
        y_out(n, s) = y_val;
        
        % 累加步长
        if cfg.useNLMS
            avg_power = local_power / max(1, (nr * ne));
            total_mu = total_mu + cfg.muMax / (avg_power + 1e-12);
        else
            total_mu = total_mu + cfg.muMax;
        end
    end
end

% === Step 4: 调试信息归一化 ===
dbg.weightNorm = sqrt(total_weight_norm_sq / max(1, (N * nr * ne)));  % 每扬声器平均 L2 范数
dbg.muUsed = total_mu / max(1, (N * ns));  % 每样本每扬声器平均步长

% === Step 5: 更新历史缓冲（用于下一帧）===
if Lw > 1
    for s = 1:ns
        for r = 1:nr
            for e_idx = 1:ne
                if N >= Lw
                    st.xfHistory(:, r, s, e_idx) = Xf(end-Lw+2:end, r, s, e_idx);
                else
                    old_hist = st.xfHistory(:, r, s, e_idx);
                    new_hist = [Xf(:, r, s, e_idx); old_hist];
                    st.xfHistory(:, r, s, e_idx) = new_hist(1:Lw-1);
                end
            end
        end
    end
end

end

% === 辅助函数：安全获取 cfg 字段 ===
function val = getfieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end