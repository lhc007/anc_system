function st = fxlms_recursive_init(cfg, S, delayD, nr, ne, ns)
% 优化版 FXLMS 初始化（保持原接口 + 支持多误差通道 + 安全配置）
%
% 输入：
%   cfg         : 配置结构体
%   S           : [Ls_max × ne × ns] 次级路径估计（允许各通道不同长度）
%   delayD      : 延迟样本数（用于参考信号对齐）
%   nr/ne/ns    : 通道数
%
% 输出：
%   st          : 初始化状态结构体（与 anc_main.m 兼容）

L = cfg.timeFilterLen;  % 控制滤波器长度

% === Step 1: 权重 W —— 支持多误差通道 (r,s,e) ===
% 使用 cell{r,s,e} 以支持独立更新（与 process 一致）
st.W = cell(nr, ns, ne);
for e = 1:ne
    for s = 1:ns
        for r = 1:nr
            st.W{r, s, e} = zeros(L, 1);
        end
    end
end

% === Step 2: 次级路径 S —— 转为 cell{e,s} 支持变长 IR ===
st.S = cell(ne, ns);
for e = 1:ne
    for s = 1:ns
        ir_raw = squeeze(S(:, e, s));  % [Ls × 1]
        % 移除尾部零（可选，提升 filter 效率）
        lastNonZero = find(ir_raw ~= 0, 1, 'last');
        if ~isempty(lastNonZero)
            ir_raw = ir_raw(1:lastNonZero);
        else
            ir_raw = ir_raw(1:1); % 全零则保留一个样本
        end
        st.S{e, s} = ir_raw(:);  % 确保列向量
    end
end

% === Step 3: 延迟与参考缓冲 ===
st.delayD = delayD;
bufferLen = delayD + L;
st.refBuffer = zeros(bufferLen, nr);  % 环形缓冲
st.refIdx = 1;                        % 下一写入位置（1-based）

% === Step 4: Filtered-X 缓冲（保留扩展性，但默认禁用）===
% 若 process 使用 conv/filter 批量计算，则 xfRing 可不用
st.xfRing = [];   % 显式清空，节省内存
st.xfIdx  = [];   % 与 xfRing 一致

% === Step 5: FIR 滤波器初始状态 zi{r,s,e}（用于 filter 模式）===
st.zi = cell(nr, ns, ne);
for r = 1:nr
    for s = 1:ns
        for e = 1:ne
            Ls = length(st.S{e, s});
            st.zi{r, s, e} = zeros(max(0, Ls - 1), 1);
        end
    end
end

% === Step X: 初始化 xfHistory（用于跨帧 Filtered-X 历史）===
if cfg.timeFilterLen > 1
    Lw = cfg.timeFilterLen;
    st.xfHistory = zeros(Lw-1, nr, cfg.numSpeakers, ne); % [Lw-1 x nr x ns x ne]
else
    st.xfHistory = [];
end

% === Step 6: 配置参数（安全读取）===
st.mu = getfieldOrDefault(cfg, 'muInit', 0.01);
st.useNLMS = getfieldOrDefault(cfg, 'useNLMS', false);
st.nlmsEps = 1e-8;
st.weight_decay = getfieldOrDefault(cfg, 'weight_decay', 0);
st.maxWeightAbs = getfieldOrDefault(cfg, 'maxWeightAbs', inf);

% === Step 7: 自适应控制 ===
st.adaptEnable = false;

% === Step 8: 反馈路径状态（由 anc_main.m 初始化）===
st.zi_fb = [];  % 占位，外部初始化

end

% === 辅助函数：安全获取字段 ===
function val = getfieldOrDefault(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end