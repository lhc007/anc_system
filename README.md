# 管道 ANC（时域递推卷积 + Xf环形缓冲）使用说明

## 文件
- anc_config.m：配置生成器（支持任意字段覆盖）
- check_anc_dependencies.m：依赖检查（输入噪声、路径文件、函数路径）
- record_secondary_path.m：多次扫频录制次级路径（对齐、异常剔除、截断）
- fxlms_recursive.m：核心算法（初始化 + 递推处理）
- anc_plot_results.m：绘图工具
- anc_main.m：主仿真入口（延迟 FXLMS）

## 流程
1. 根据设备与文件名修改 anc_config.m。
2. 运行 record_secondary_path.m 生成 secondary_path.mat。
3. 准备 feedback_path.mat（变量 F），尺寸 [Lfb x numRef x numSpeakers]。
4. 准备多通道输入噪声 anc_sim_data_road_noise.wav（与 cfg.outputFile 匹配）。
5. 运行 anc_main.m 启动仿真。

## 算法要点
- 延迟 D = 次级路径峰值中位数 + 安全裕量（cfg.delayMarginSamples）。
- 每个 (r,s,e) 建立一个 FIR 递推滤波状态（filter + zi），每新样本输出一个 Filtered-X 样本。
- 维护一个长度 L 的环形缓冲（xfRing）保存最近 L 个 Filtered-X 样本，更新权重时直接读取，避免重复卷积。
- 梯度 = Σ_e e_e(n) * Xf_e（逐误差通道累加）。
- 可选 NLMS（useNLMS=true），默认使用固定步长调度（更易对比与控制）。

## 调优
- 提升收敛：增大 muMax，减小 L；合理延迟 D（避免过大导致预测滞后）。
- 稳定性：增大 weight_decay，减小 maxWeightAbs。
- 目标低频：在录制中关闭 secUseFreqWeight 或调低 lowBoostHz。
- 异常：检查 S 的有效能量是否在前 512~1024 样本内，必要时调小 irTruncateLen。

## 注意
- 本方案优先服务于短管道、短时延场景，效率高、延迟小。若 L 很大（>2048）或 fs 更高，可考虑频域 Block FXLMS 或混合方案。
- 若需要反馈路径录制脚本，请告知，我会提供 record_feedback_path.m。

祝实验顺利！