function cfg = anc_config(varargin)
% ANC 配置生成器（管道主动降噪系统） 支持任意字段用 Name, Value 形式覆盖:
%   cfg = anc_config('fs', 48000, 'timeFilterLen', 1024);
%
% 系统描述：
%   管道长度约 40 cm，布置顺序：参考麦阵列 (4) -> 扬声器阵列 (4) -> 误差麦阵列 (2)
%   两端开口，主要控制低频路噪（胎噪、发动机低阶、风噪低频部分）。
%
% 重要字段说明：
%   fs                 采样率 frameSize          每帧处理样本数 timeFilterLen
%   控制滤波器长度（FXLMS 自适应 FIR 长度） secondaryPathFile  次级路径文件（录制脚本生成）
%   feedbackPathFile   反馈路径文件（扬声器→参考麦泄漏路径） outputFile
%   离线仿真输入噪声的多通道文件 maxOutput          扬声器输出限幅
%
% 模式：
%   useNLMS            是否使用 NLMS（归一化步长），默认 false enableMinPhaseInit
%   是否用最小相位次级路径初始化权重 enableSecPathWeight 是否在录制时使用低频权重
%
% 日志：
%   logging            是否保存运行日志 logFile            日志保存文件名
%
% 自适应时序：
%   startMuteFrames    初始静音帧（不输出、不适应） rampFrames         输出渐增帧 warmupFrames
%   步长升温帧（从初始值到最大值）
%
% 步长：
%   muInit / muMax / muMin    步长范围（若 useNLMS=false 则固定调度）
%
% 延时：
%   delayMarginSamples        基于次级路径峰值自动估计延时时添加的安全裕量
%
% 能量保护：
%   weight_decay              权重衰减系数（正则） maxWeightAbs              权重幅值限制
%
% 监控带通：
%   bandFreqLow / bandFreqHigh 用于监控改善量的频段（不影响自适应）
%
% 作者：刘昊诚

cfg = struct();

%% 基础系统参数
cfg.fs               = 48000;       %16000->48000
cfg.frameSize        = 128;         % 小帧以降低延迟；若收敛慢可增大到 128或256
cfg.timeFilterLen    = 3072;         % 控制滤波器长度；若次级路径较长，建议 128 或 256
cfg.numSpeakers      = 2;       %先用两个扬声器做测试效果合格在改为四个
% 限制 ANC 带宽（聚焦马路噪音）
cfg.ancLowpassHz = 1000;   % 新增字段：ANC 工作上限频率

% 硬件设备
cfg.micDeviceName          = '六通道麦克风阵列 (YDM6MIC Audio)';
cfg.spkDevice1Name         = '扬声器 (USB Audio Device)';
cfg.spkDevice2Name         = '扬声器2 (Realtek(R) Audio)';
% 扬声器映射简化：只指定通道号（1~4）
cfg.spkChannels = [1, 2, 3, 4];  % Spk1→Ch1, Spk2→Ch2, Spk3→Ch3, Spk4→Ch4

%% 扬声器映射（关键！）
% 每个扬声器指定其输出设备和通道 格式: {设备名称, 输出通道索引 (1=左, 2=右)}
cfg.spkMapping = {
    {'扬声器 (USB Audio Device)', 1};   % 扬声器 1 → USB 左(左上)
    {'扬声器 (USB Audio Device)', 2};   % 扬声器 2 → USB 右(左下）
    {'扬声器2 (Realtek(R) Audio)', 1};  % 扬声器 3 → Realtek 左(右上)
    {'扬声器2 (Realtek(R) Audio)', 2}   % 扬声器 4 → Realtek 右(右下)
    };


%% 麦克风与通道映射
cfg.micNumChannels         = 6;
cfg.micChannels.reference  = [1 2 3 4];
cfg.micChannels.error      = [5 6];

%% 文件路径（明确区分输入/输出）
cfg.inputAudioFile    = 'noise/anc_sim_data_road_noise.wav';        % ← 输入噪声（多通道）
cfg.secondaryPathFile = 'secondary_path/secondary_path.mat';           % 次级路径地址
cfg.feedbackPathFile  = 'feedback_path/feedback_path.mat';            % 反馈路径
cfg.logFile           = 'anc_run_log.mat';              %运行日志
cfg.outputAudioFile   = 'anc_output.wav';   % ← 仿真输出
cfg.noiseFile = 'road_noise.wav';           %初始噪音文件


%% 次级路径录制（ 测量专用参数）

cfg.spkAmplitude = [1.2, 1.2]; % 默认为 1.0

cfg.repetitions = 4;               % 增加重复测量次数，提高统计可靠性
cfg.timeFrameSamples = 512;         % 帧大小 必须是 2 的幂，且与 ASIO buffer 匹配
cfg.irMaxLen              = 4096;
cfg.preRollFrames         = 20;  % 预热帧
cfg.tailNoiseLen          = 512;
cfg.saveFirstRaw          = true;
cfg.saveAllRaw            = true;    %保存次级路径原始录音
cfg.doAlignRepeats        = false;   % 保留真实延迟
cfg.exportAlignedIR       = true;    % 诊断用对齐 IR
cfg.preSilenceSec         = 0.8;
cfg.postSilenceSec        = 0.4;
cfg.writeBlockPad         = true;

% 峰与漂移相关参数
cfg.enableLowFreqBoost       = true;
cfg.lowFreqCutHz             = 150;
cfg.lowFreqMixRatio          = 0.9;     % 90% 能量集中在低频
cfg.minPhysDelaySamples      = 50;
cfg.maxAllowedDriftSamples   = 200;
cfg.enableRepeatAlignment    = true;
% SNR 排除：尝试剔除低 SNR 重复
cfg.excludeLowSNR            = true;
cfg.snrThresholdDB           = 6;

% 可靠峰判定参数
cfg.reliableMinRatio        = 0.30;         % 可靠比例阈值
cfg.reliableMaxIQR          = 1200;          % 适当放宽IQR阈值
cfg.filterTailFraction      = 0.5;
cfg.enableOutlierReject     = true;
cfg.outlierIQRMultiplier    = 1.2;    % 更严格地剔除离群点

% deconvolve_sweep 参数（可选暴露）
cfg.deconvExtraTail        = cfg.irMaxLen + 3000;
cfg.prePeakKeep            = 32;
cfg.deconvPreDelayKeep     = 768;
cfg.deconvPeakThreshDB     = 3;      % 从 5 -> 3（低 SNR 下更易检测峰值）
cfg.deconvMaxSearch        = 9000;
cfg.deconvRegEps           = 1e-3;  % 从 1e-4 -> 1e-3 提升数值稳定性（避免放大噪声）
cfg.deconvNoiseWin         = 1200;  % 从 600 -> 1200，长窗更好估计噪声
cfg.envSmoothWin           = 10;
cfg.deconvCumEnergyFrac    = 0.05;
cfg.deconvMinPeakFrac      = 0.005;
cfg.deconvSnrBodyRadius    = 96;
cfg.deconvFftCorrEnable    = true;

% sweepCfg参数
cfg.padLeading = 0.25;
cfg.padTrailing = 0.25;
cfg.amplitude = 0.98;
cfg.sweepDuration = 16;

%% ========== 反馈路径测量参数 (Feedback Path) ==========
% 扫频信号设置
cfg.fbSweepDur        = 3.2;          % 扫频时长 (秒) → 对应 153600 样本 @48kHz
cfg.fbSweepFstart     = 100;          % 起始频率 (Hz)，避开极低频噪声
cfg.fbSweepFend       = 8000;         % 截止频率 (Hz)，覆盖主要控制带宽
cfg.fbSweepNsweeps    = 2;            % 重复测量次数，用于平均降噪
cfg.feedbackIRLenSec = 0.5;           % 500 毫秒

% 频率加权（增强低频信噪比）
cfg.fbUseFreqWeight   = true;         % 启用低频提升
cfg.fbLowBoostHz      = 250;          % 低于此频率开始加权
cfg.fbLowBoostPower   = 0.6;          % 加权斜率 (≈ 1/sqrt(f))，避免过载

% 峰值对齐与可靠性
cfg.fbAlignTargetIdx  = [];           % 自动选择主峰（留空）
cfg.fbAlignOffset     = 12;           % 主峰前保留样本数（因果保护）
cfg.fbAlignThreshDb   = -30;          % 峰值检测阈值（dB，相对于峰值）
cfg.fbCorrMin         = 0.75;         % 重复间互相关下限，剔除异常

% 冲激响应截断
cfg.fbEnergyCutRatio  = 0.9995;       % 保留 99.95% 能量（平衡长度与精度）
cfg.irMaxLen          = 4096;         % 最大 IR 长度（样本），@48kHz ≈ 85ms
cfg.irTruncateLen     = 0;            % 0 表示自动截断；>0 则强制固定长度

%% 时序与步长
cfg.startMuteFrames        = 20;
cfg.rampFrames             = 300;
cfg.warmupFrames           = 200;

%% 步长与算法选项
cfg.muInit                 = 1e-6;
cfg.muMax                  = 5e-5;  % 若系统稳定，可加速收敛1e-4
cfg.muMin                  = 1e-7;

cfg.useNLMS                = false;
cfg.enableMinPhaseInit     = false;

%% 延迟与稳定性
cfg.delayMarginSamples     = 16;      % 若出现相位/延迟抖动，可提升到 8~16
cfg.weight_decay           = 1e-5;
cfg.maxWeightAbs           = 0.05;

%% 反馈泄漏相关（是否使用反馈泄漏抵消功能）
cfg.useFeedbackLeakCancel  = false;

%% 监控频段
cfg.bandFreqLow            = 20;
cfg.bandFreqHigh           = 800;

%% 输出限幅与保护
cfg.maxOutput              = 0.3;

%% 运行日志
cfg.statusPrintFrames      = 100;
cfg.logging                = true;

%% 离线仿真相关
cfg.duration_sec           = 30;

%% 噪声基线估计
cfg.noiseBaselineFrames    = 80;   % 至少 80 帧观察
cfg.noiseBaselineMinCount  = 40;   % 至少 40 个样本才锁定

%% 在线建模与自适应对齐
cfg.enableOnlineSecPath    = true;
cfg.probeAmplitude         = 1e-4; % 探测信号幅度；约 -50 dBFS 相对 0.3 限幅
cfg.enableAdaptiveDelay    = true;
cfg.minDelaySamples        = 8;
cfg.maxDelaySamples        = 32;
cfg.delaySearchRange       = [64, 256];
cfg.enableBeamforming      = true;
cfg.refArraySpacing_m      = 0.02;

%% 可选字段覆盖
if mod(nargin,2)~=0
    error('anc_config: 必须成对传入参数 (Name,Value)');
end
for k = 1:2:nargin
    name = varargin{k};
    if isfield(cfg, name)
        cfg.(name) = varargin{k+1};
    else
        warning('anc_config: 未知配置字段 "%s" 被忽略', name);
    end
end
end