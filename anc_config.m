function cfg = anc_config(varargin)
% ANC 配置生成器（管道主动降噪系统） 支持任意字段用 Name, Value 形式覆盖:
%   cfg = anc_config('fs', 48000, 'timeFilterLen', 1024);
%
% 系统描述：
%   管道长度约 50 cm，布置顺序：参考麦阵列 (4) -> 扬声器阵列 (4) -> 误差麦阵列 (2)
%   两端开口，主要控制低频路噪（胎噪、发动机低阶、风噪低频部分）。
%
% 重要字段说明：
%   fs                 采样率 
%   frameSize          每帧处理样本数 
%   timeFilterLen 控制滤波器长度（FXLMS 自适应 FIR 长度） 
%   secondaryPathFile  次级路径文件（录制脚本生成）
%   feedbackPathFile   反馈路径文件（扬声器→参考麦泄漏路径） 
%   inputAudioFile 离线仿真输入噪声的多通道文件 
%   maxOutput          扬声器输出限幅
%
% 模式：
%   useNLMS            是否使用 NLMS（归一化步长），默认 false 
%   enableMinPhaseInit 是否用最小相位次级路径初始化权重 enableSecPathWeight 是否在录制时使用低频权重
%
% 日志：
%   logging            是否保存运行日志 logFile            日志保存文件名
%
% 自适应时序：
%   startMuteFrames    初始静音帧（不输出、不适应） 
%   rampFrames         输出渐增帧 
%   warmupFrames 步长升温帧（从初始值到最大值）
%
% 步长：
%   muInit / muMax / muMin    步长范围（若 useNLMS=false 则固定调度）
%
% 延时：
%   delayMarginSamples        基于次级路径峰值自动估计延时时添加的安全裕量
%
% 能量保护：
%   weight_decay              权重衰减系数（正则） 
%   maxWeightAbs              权重幅值限制
%
% 监控带通：
%   bandFreqLow / bandFreqHigh 用于监控改善量的频段（不影响自适应）
%
% 作者：刘昊诚

cfg = struct();

%% 基础系统参数
cfg.fs               = 48000;       % 采样率（Hz），决定系统最高可处理频率（奈奎斯特频率为 24 kHz）
cfg.frameSize        = 128;         % 每帧处理的样本数；影响算法延迟和实时性（128 @48kHz ≈ 2.7 ms/帧）
cfg.timeFilterLen    = 3072;        % FXLMS 自适应滤波器的抽头数（FIR 长度）；需覆盖次级路径主能量持续时间
cfg.numSpeakers      = 2;           % 当前激活的扬声器数量（最大支持 4，先用 2 个测试）
cfg.ancLowpassHz     = 1000;        % ANC 系统工作上限频率（高于此频段不进行降噪处理）

%% ANC系统监控频段
cfg.bandFreqLow      = 60;          % 用于计算降噪效果（如 dB 改善量）的频带下限（Hz）
cfg.bandFreqHigh     = 800;         % 用于计算降噪效果的频带上限（Hz）；覆盖典型路噪主能量区

%% 输出限幅与保护
cfg.maxOutput        = 0.3;         % 扬声器输出信号的最大幅值（防止削波或硬件过载，0.3 ≈ -10.5 dBFS）

% 硬件设备
cfg.micDeviceName    = '六通道麦克风阵列 (YDM6MIC Audio)';  % 多通道麦克风设备名称（Windows 音频设备）
cfg.spkDevice1Name   = '扬声器 (USB Audio Device)';         % 第一个扬声器输出设备（USB 声卡）
cfg.spkDevice2Name   = '扬声器2 (Realtek(R) Audio)';        % 第二个扬声器输出设备（主板 Realtek 声卡）

%% 扬声器映射（关键！）
% 每个扬声器指定其输出设备和通道 格式: {设备名称, 输出通道索引 (1=左, 2=右)}
% cfg.spkMapping = {
%     {'扬声器 (USB Audio Device)', 1};   % 扬声器 1 → USB 声卡左声道（通常对应物理左上位置）
%     {'扬声器 (USB Audio Device)', 2};   % 扬声器 2 → USB 声卡右声道（通常对应物理左下位置）
%     {'扬声器2 (Realtek(R) Audio)', 1};  % 扬声器 3 → Realtek 左声道（通常对应物理右上位置）
%     {'扬声器2 (Realtek(R) Audio)', 2}   % 扬声器 4 → Realtek 右声道（通常对应物理右下位置）
% };

%% 麦克风与通道映射
cfg.micNumChannels         = 6;                % 麦克风阵列总通道数
cfg.micChannels.reference  = [1 2 3 4];        % 参考麦克风通道索引（用于获取初级噪声参考信号）
cfg.micChannels.error      = [5 6];            % 误差麦克风通道索引（用于反馈控制误差信号）
cfg.numRefs = 4;                               % 参考麦克风数量
cfg.numErrs = 2;                               % 误差克风数量

%% 文件路径（明确区分输入/输出）
cfg.inputAudioFile    = 'noise/anc_sim_data_road_noise.wav';        % 输入的多通道路噪仿真信号（用于离线测试）
cfg.secondaryPathFile = 'secondary_path/secondary_path.mat';        % 次级路径冲激响应文件（由 record_secondary_path.m 生成）
cfg.feedbackPathFile  = 'feedback_path/feedback_path.mat';         % 反馈路径（扬声器→参考麦）冲激响应文件
cfg.logFile           = 'anc_run_log.mat';                         % 运行过程日志（含权重、SNR、改善量等）
cfg.outputAudioFile   = 'anc_output.wav';                          % ANC 处理后的输出音频（多通道或单通道）
cfg.noiseFile         = 'road_noise.wav';                          % 原始单通道路噪样本（用于生成仿真数据）

%% 次级路径录制（测量专用参数）

cfg.saveEachRepeatIR = true;      % 保存每个repeat的IR用于后期分析
cfg.saveDiagnosticInfo = true;     % 保存完整的元数据
cfg.enableAutoGainAdjustment =true; %SNR过低时给出增益调整建议

% sweepCfg参数（扫频信号生成）
cfg.padLeading        = 0.5;    % 扫频前导静音时间（秒）
cfg.padTrailing       = 1;      % 扫频尾随静音时间（秒）
cfg.amplitude         = 0.98;   % 扫频信号幅值（接近满幅但避免削波）
cfg.sweepDuration     = 5;      % 扫频持续时间（秒）；越长频率分辨率越高，低频能量越强
cfg.minSnrForReliable = 3; 
cfg.spkAmplitude      = [0.9, 0.9];            % 播放扫频信号时各扬声器的增益（>1 表示数字域放大，需注意不削波）
cfg.sweepF1           = 80;                    % 扫频信号起始频率（Hz）；避开无效低频（<60 Hz）
cfg.sweepF2           = 1200;                  % 扫频信号终止频率（Hz）；覆盖 ANC 主要工作带宽
cfg.repetitions       = 1;                     % 重复播放并录制次数；用于提高 SNR 和鲁棒性
cfg.timeFrameSamples  = 1024;                  % 录音/播放缓冲区帧大小（必须是 2 的幂）
cfg.irMaxLen          = 4096;                 % 冲激响应最大截断长度（样本数，@48kHz ≈ 85 ms）
cfg.preRollFrames     = 20;                    % 开始正式记录前的预热帧数（丢弃初始不稳定数据）
cfg.tailNoiseLen      = 512;                   % 用于估计噪声底噪的尾部静音段长度（样本）
cfg.saveFirstRaw      = true;                 % 是否保存第一次原始录音（用于调试）
cfg.saveAllRaw        = true;                  % 是否保存所有重复的原始录音（占用磁盘空间）
cfg.doAlignRepeats    = true;                  % 是否对多次重复测量进行时间对齐（false 保留真实延迟变化）
cfg.exportAlignedIR   = true;                  % 是否导出对齐后的平均 IR（用于诊断）
cfg.preSilenceSec     = 0.8;                   % 扫频前静音时间（秒），确保系统稳定
cfg.postSilenceSec    = 0.4;                   % 扫频后静音时间（秒），捕获完整 IR 尾部
cfg.writeBlockPad     = false;                 % 是否在播放块前后填充静音（避免突发）

% 峰与漂移相关参数
cfg.enableLowFreqBoost    = false;             % 是否在扫频中增强低频能量（提升低频 SNR）
cfg.lowFreqCutHz          = 80;                % 低频增强的截止频率（低于此频率按比例提升）
cfg.lowFreqMixRatio       = 0.8;               % 低频增强混合比例（0.9 表示 90% 能量集中在低频段）
cfg.maxAllowedDriftSamples = 200;              % 允许的最大重复间延迟漂移（样本）；超限则标记不可靠
cfg.enableRepeatAlignment = true;              % 是否启用重复测量间的自动对齐（基于互相关）
cfg.enableRealTimeMonitor = true;
% SNR 排除：尝试剔除低 SNR 重复
cfg.excludeLowSNR         = true;              % 是否排除 SNR 低于阈值的重复测量
% SNR阈值
cfg.snrThresholdDB = 10;                       % SNR 可靠性阈值（dB）；低于此值视为不可用

% 可靠峰判定参数
cfg.reliableMinRatio      = 0.30;              % 可靠重复占比阈值（如 4 次中有 ≥2 次可靠才接受）
cfg.reliableMaxIQR        = 1200;              % 延迟估计的四分位距（IQR）上限（样本）；过大表示分散
cfg.filterTailFraction    = 0.5;               % 在尾部加窗时保留的比例（用于平滑 IR 截断）
cfg.enableOutlierReject   = true;              % 是否启用离群点剔除（基于 IQR）
cfg.outlierIQRMultiplier  = 1.2;               % IQR 倍数阈值；超出 [Q1 - k*IQR, Q3 + k*IQR] 视为离群

% deconvolve_sweep 参数（反卷积扫频恢复 IR）
cfg.deconvExtraTail       = cfg.irMaxLen + 3000; % 反卷积时额外补零长度（提升频率分辨率）
cfg.prePeakKeep           = 256;                  % 主峰前保留的样本数（用于因果性检查）
cfg.deconvPreDelayKeep    = 768;                 % 反卷积结果中强制保留的前导静音长度
cfg.deconvMaxSearch       = 20000;                % 最大搜索延迟范围（样本）
cfg.deconvRegEps          = 1e-4;                % Tikhonov 正则化参数；抑制噪声放大（值越大越平滑）
cfg.deconvNoiseWin        = 800;                % 用于估计噪声功率的窗口长度（样本）
cfg.deconvEnvSmoothWin    = 10;                  % 包络平滑窗口长度（样本）
cfg.deconvSnrBodyRadius   = 96;                  % SNR 计算时主峰邻域半径（样本）
cfg.deconvFftCorrEnable   = true;                % 是否启用 FFT 域互相关辅助峰值定位
cfg.deconvDebugMode       = true;                % ✅ 新增：启用deconvolve_sweep的调试模式

% 在 anc_config.m 中添加
cfg.minPhysDelaySamples = round(0.002 * cfg.fs); % 2ms = 96样本 最小物理延迟（样本）
cfg.maxPhysDelaySamples = round(00.1 * cfg.fs);    % 100ms = 4800样本 最大物理延迟（样本）
cfg.delaySearchRadius = 2000;   % 延迟搜索半径（样本）
cfg.peakRefineRadius = 150;     % 峰值细化半径（样本）

% 预回声容忍度分级
cfg.preEchoSevereThresh = 0.15; % 严重预回声阈值
cfg.preEchoModerateThresh = 0.05; % 中等预回声阈值

% 增加以下参数
cfg.enableIRDiagnostic = true;  % 启用IR诊断数据保存
cfg.peakDetectionMinDistance = 100;  % 峰值最小距离
cfg.peakDetectionProminence = 0.1;   % 峰值显著性阈值

% 修改反卷积参数
cfg.deconvPeakThreshDB = 12;     % ✅ 提高峰值检测阈值，从12提高到15
cfg.deconvCumEnergyFrac = 0.05;  % ✅ 提高累积能量阈值，从0.05提高到0.10
cfg.deconvMinPeakFrac = 0.005;   % ✅ 提高最小峰值比例，从0.005提高到0.008

%% ========== 反馈路径测量参数 (Feedback Path) ==========
% 扫频信号设置
cfg.fbSweepDur        = 3.2;          % 反馈路径扫频时长（秒）
cfg.fbSweepFstart     = 100;          % 反馈路径扫频起始频率（Hz），避开极低频噪声
cfg.fbSweepFend       = 8000;         % 反馈路径扫频终止频率（Hz），覆盖控制带宽
cfg.fbSweepNsweeps    = 2;            % 反馈路径重复测量次数
cfg.feedbackIRLenSec  = 0.5;          % 反馈路径 IR 最大长度（秒）

% 频率加权（增强低频信噪比）
cfg.fbUseFreqWeight   = true;         % 是否对反馈扫频进行低频加权
cfg.fbLowBoostHz      = 250;          % 低频加权起始频率（Hz）
cfg.fbLowBoostPower   = 0.6;          % 加权斜率（≈ f^{-0.6}），提升低频 SNR 但避免过载

% 峰值对齐与可靠性
cfg.fbAlignTargetIdx  = [];           % 指定对齐目标通道（空表示自动选主峰）
cfg.fbAlignOffset     = 12;           % 对齐后主峰前保留的样本数（保证因果性）
cfg.fbAlignThreshDb   = -30;          % 峰值检测阈值（dB）
cfg.fbCorrMin         = 0.75;         % 重复测量间互相关系数下限；低于则剔除

% 冲激响应截断
cfg.fbEnergyCutRatio  = 0.9995;       % 保留 IR 总能量的 99.95% 以确定截断点
cfg.irTruncateLen     = 0;            % 强制截断长度（0 表示自动）

%% 时序与步长
cfg.startMuteFrames   = 20;           % 启动后静音帧数（不输出、不更新权重）
cfg.rampFrames        = 300;          % 输出幅值从 0 线性上升到 maxOutput 的帧数（防冲击）
cfg.warmupFrames      = 200;          % 步长从 muInit 升至 muMax 的过渡帧数

%% 步长与算法选项
cfg.muInit            = 1e-6;         % 初始步长（FXLMS 算法学习率）
cfg.muMax             = 5e-5;         % 最大步长（收敛快但可能不稳定）
cfg.muMin             = 1e-7;         % 最小步长（防发散）
cfg.useNLMS           = false;        % 是否使用归一化 LMS（NLMS）；true 时步长自动归一化
cfg.enableMinPhaseInit= false;        % 是否用最小相位近似初始化控制器权重（加速收敛）

%% 延迟与稳定性
cfg.delayMarginSamples= 16;           % 在估计的次级路径延迟基础上增加的安全裕量（样本）
cfg.weight_decay      = 1e-5;         % 权重衰减系数（L2 正则化，防过拟合和漂移）
cfg.maxWeightAbs      = 0.05;         % 权重绝对值上限（防数值爆炸）

%% 反馈泄漏相关（是否使用反馈泄漏抵消功能）
cfg.useFeedbackLeakCancel = false;    % 是否启用反馈泄漏抵消（当参考麦收到扬声器泄漏时启用）

%% 运行日志
cfg.statusPrintFrames = 100;          % 每隔多少帧打印一次状态信息（如改善量）
cfg.logging           = true;         % 是否保存详细运行日志到 cfg.logFile

%% 离线仿真相关
cfg.duration_sec      = 30;           % 离线仿真运行时长（秒）

%% 噪声基线估计
cfg.noiseBaselineFrames   = 80;       % 用于估计背景噪声能量的初始帧数
cfg.noiseBaselineMinCount = 40;       % 至少需要多少帧有效噪声样本才锁定基线

%% 在线建模与自适应对齐
cfg.enableOnlineSecPath   = true;     % 是否在线估计次级路径（应对环境变化）
cfg.probeAmplitude        = 1e-4;     % 在线探测信号幅值（叠加在输出上，用于辨识）
cfg.minDelaySamples       = 8;        % 允许的最小延迟（样本）
cfg.maxDelaySamples       = 32;       % 允许的最大延迟（样本）
cfg.delaySearchRange      = [64, 256];% 在线延迟搜索范围（样本）
cfg.enableBeamforming     = true;     % 是否对参考麦克风阵列做波束成形（聚焦前方噪声源）
cfg.refArraySpacing_m     = 0.02;     % 参考麦克风阵元间距（米），用于波束成形计算

% === 自适应延迟（Adaptive Delay）相关配置 ===
% 启用/禁用自适应延迟模块
cfg.enableAdaptiveDelay = false; 

% 初始等待与触发控制（以帧为单位）
cfg.adaptiveDelayWarmupFrames = 500; 
cfg.adaptiveDelayTriggerFrames = 200; 
% 每次调用 adaptive_delay_estimator 内部是否做估计的帧周期（estimator 内部也有 update 控制）
cfg.adaptiveDelayUpdateFrames = 50; 

% GCC-PHAT 相关：用于互相关计算的窗口长度（秒）
cfg.adaptiveDelayWindowSec = 0.2;

% 置信度门限（峰值/旁瓣比），越大越严格
cfg.adaptiveDelayConfidenceThreshold = 6.0; 

% EMA 平滑因子（0~1），越大响应越快但更抖动
cfg.adaptiveDelaySmoothAlpha = 0.3; 

% 单次允许的最大/最小延迟变化（样本）
cfg.adaptiveDelayMaxChangeSamples = round(0.02 * cfg.fs); 
cfg.adaptiveDelayMinChangeSamples = 1; 

% 初始延迟（可覆盖来自次级路径文件的估计）
cfg.initialDelaySamples = round(0.0012 * cfg.fs); 

% 置信度/更新行为调优（可选）
% 是否只在 probe/高 SNR 段允许更新（默认 false，可在主循环中按需检测）
cfg.adaptiveDelayRequireHighSNR = false; 
cfg.adaptiveDelayMinRefRms = 1e-6; 
cfg.adaptiveDelayMinErrRms = 1e-6; 

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