function st = fxlms_recursive_init(cfg, S, delayD, nr, ne, ns)
% 优化版 FXLMS 初始化（数组化 W + 次级路径 cell 化 + 有状态滤波支持）
%
% 输入：
%   cfg         : 配置结构体（含 timeFilterLen、muInit 等）
%   S           : [Ls_max × ne × ns] 次级路径（允许各通道不同长度，尾零可存在）
%   delayD      : 延迟样本数（用于参考对齐）
%   nr/ne/ns    : 通道数
%
% 输出：
%   st          : 初始化状态结构体

Lw = cfg.timeFilterLen;

% ---------------- 权重矩阵（数组）----------------
% W: [Lw × nr × ns × ne]，全零初始化
st.W = zeros(Lw, nr, ns, ne);

% ---------------- 次级路径（cell 化）----------------
% st.S_cell{e,s} = 列向量 IR（移除尾零）
st.S_cell = cell(ne, ns);
for e = 1:ne
    for s = 1:ns
        ir_raw = squeeze(S(:, e, s));  % [Ls × 1]
        if ~isvector(ir_raw)
            ir_raw = ir_raw(:);
        end
        lastNonZero = find(ir_raw ~= 0, 1, 'last');
        if ~isempty(lastNonZero)
            ir = ir_raw(1:lastNonZero);
        else
            ir = 0; % 全零时保留一个样本
        end
        st.S_cell{e, s} = ir(:);
    end
end

% ---------------- 参考环形缓冲 ----------------
st.delayD   = delayD;
st.refBuffer = zeros(delayD + Lw, nr);
st.refIdx    = 1;  % 下一个写入位置（1-based）

% ---------------- Filtered-X 历史窗口 ----------------
% 用于构造长度 Lw 的滑动窗口（向量化）
st.xfHistory = zeros(max(Lw - 1, 0), nr, ns, ne);

% ---------------- filter+zi 的状态（可选）----------------
% 当 fxlmsRecursiveProcess 使用 'filter' 模式时需要：
% st.zi_fxf{e,s,r} 长度为 length(st.S_cell{e,s}) - 1
st.zi_fxf = cell(ne, ns, nr);
for e = 1:ne
    for s = 1:ns
        Ls = length(st.S_cell{e, s});
        for r = 1:nr
            st.zi_fxf{e, s, r} = zeros(max(0, Ls - 1), 1);
        end
    end
end

% ---------------- 反馈路径状态（外部初始化）----------------
st.zi_fb = []; % 在主程序中初始化为 cell(numRef, numSpeakers)

% ---------------- 自适应控制与配置 ----------------
st.adaptEnable  = false;
st.mu           = getfieldOrDefault(cfg, 'muInit', 0.01);
st.useNLMS      = getfieldOrDefault(cfg, 'useNLMS', false);
st.nlmsEps      = 1e-8;
st.weight_decay = getfieldOrDefault(cfg, 'weight_decay', 0);
st.maxWeightAbs = getfieldOrDefault(cfg, 'maxWeightAbs', inf);

end

% 安全获取字段
function val = getfieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end