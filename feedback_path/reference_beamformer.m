function ref_bf = reference_beamformer(ref_raw, fs, mic_pos, look_dir_deg)
% 简单 delay-and-sum 波束成形
% 线性麦克风阵列 delay-and-sum 波束成形
% ref_raw: [N x nr]
% mic_pos: 麦克风位置向量（单位：米，相对原点）
% look_dir_deg: 期望方向（0° = 正对管道轴向，从室外来）

if nargin < 4, look_dir_deg = 0; end
c = 343; % 声速 m/s
theta = deg2rad(look_dir_deg);
nr = size(ref_raw,2);
if nargin < 3 || numel(mic_pos) ~= nr
    % 默认等距 4 麦：间距 2 cm
    mic_pos = (0:nr-1)*0.02;
end

% 计算 steering delays（以第一个麦为参考）
delays_samp = -mic_pos * sin(theta) / c * fs; % 负号：波从 theta 方向来
delays_samp = delays_samp - min(delays_samp); % 对齐到最早到达

% 简化：用整数延迟 + 插值（此处仅用整数）
ref_bf = zeros(size(ref_raw,1), 1);
for r = 1:nr
    d = round(delays_samp(r));
    if d == 0
        sig = ref_raw(:,r);
    elseif d > 0
        sig = [zeros(d,1); ref_raw(1:end-d,r)];
    else
        sig = [ref_raw(-d+1:end,r); zeros(-d,1)];
    end
    ref_bf = ref_bf + sig;
end
ref_bf = ref_bf / nr;
end