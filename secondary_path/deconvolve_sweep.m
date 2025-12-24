function out = deconvolve_sweep(recorded, sweepSig, fs, params)
% deconvolve_sweep (修复版 v1.7)
% 修复点：修复搜索范围计算错误，添加互相关延迟验证

if nargin < 4 || isempty(params), params = struct(); end

regEps        = getP(params,'regEps',1e-4);
extraTail     = getP(params,'extraTail',4096);
preDelayKeep  = getP(params,'preDelayKeep',256);
tailTotal     = getP(params,'tailTotal',4096);
peakThreshDB  = getP(params,'peakThreshDB',12);
maxSearch     = getP(params,'maxSearch',15000);
noiseWin      = getP(params,'noiseWin',400);
envSmoothWin  = getP(params,'envSmoothWin',8);
cumEnergyFrac = getP(params,'cumEnergyFrac',0.05);
minPeakFrac   = getP(params,'minPeakFrac',0.005);
snrBodyRadius = getP(params,'snrBodyRadius',96);
fftCorrEnable = getP(params,'fftCorrEnable',true);
debugMode     = getP(params,'debugMode',false);

% 物理延迟范围参数
minPhysDelay  = getP(params, 'minPhysDelay', 50);
maxPhysDelay  = getP(params, 'maxPhysDelay', 5000);
delaySearchRadius = getP(params, 'delaySearchRadius', 1000);

rec = recorded(:);
exc = sweepSig(:);
Nrec = length(rec);
Nexc = length(exc);

%% 互相关延迟估计
if fftCorrEnable
    NfftCorr = 2^nextpow2(Nrec + Nexc);
    RECf = fft(rec, NfftCorr);
    EXCf = fft(exc, NfftCorr);
    C = ifft(RECf .* conj(EXCf));
    C = [C(end-Nexc+2:end); C(1:Nrec)];
    lags = (-Nexc+2:Nrec)';
else
    [C,lags] = xcorr(rec, exc);
end
[~,imax] = max(C);
delayCorr = lags(imax);

% ✅ 修复：验证互相关延迟的合理性
if debugMode
    fprintf('  [DEB-xcorr] raw delayCorr=%d, Nrec=%d, Nexc=%d\n', ...
        delayCorr, Nrec, Nexc);
end

% 如果延迟异常大，可能是计算错误，进行修正
if abs(delayCorr) > min(Nrec, Nexc) * 0.8
    % 尝试使用绝对值互相关重新计算
    [~, imax2] = max(abs(C));
    delayCorr2 = lags(imax2);
    if debugMode
        fprintf('  [DEB-xcorr] 使用绝对值互相关，delayCorr2=%d\n', delayCorr2);
    end
    delayCorr = delayCorr2;
end

startIdx = delayCorr + 1;
warnEarly = startIdx < 1;
if startIdx < 1, startIdx = 1; end

segEnd = min(startIdx + Nexc + extraTail - 1, Nrec);
rec_aligned = rec(startIdx:segEnd);

%% 频域反卷积
Nfft = 2^nextpow2(length(rec_aligned) + Nexc - 1);
REC = fft(rec_aligned, Nfft);
EXC = fft(exc, Nfft);
magEXC2 = abs(EXC).^2;
Hf = (REC .* conj(EXC)) ./ (magEXC2 + regEps);
h_full = real(ifft(Hf));

maxLengthToKeep = max(maxSearch + 2048, 16384);
if length(h_full) > maxLengthToKeep
    h_full = h_full(1:maxLengthToKeep);
end

% 去局部均值
nw_pre = min(noiseWin, floor(length(h_full)/10));
if nw_pre < 10, nw_pre = min(100, length(h_full)); end
noiseBaseMean = mean(h_full(1:nw_pre));
h_full = h_full - noiseBaseMean;

%% 噪声统计 + 阈值
nw = min(noiseWin, length(h_full)-10);
if nw < 10, nw = min(100, length(h_full)); end
noiseSeg = h_full(1:nw);
noiseStd = std(noiseSeg);
noiseMAD = median(abs(noiseSeg - median(noiseSeg))) / 0.6745;
noiseBase = max(noiseStd, noiseMAD);
th = noiseBase * 10^(peakThreshDB/20);

searchEnd = min(length(h_full), maxSearch);
if searchEnd < 10, searchEnd = min(length(h_full), 100); end

%% ✅ 修复：正确的搜索范围计算
% 1. 首先计算包络
cand = h_full(1:searchEnd);
envRaw = abs(cand);
envSm = movmean(envRaw, envSmoothWin);
totalE = sum(envSm);
if totalE < 1e-18
    out = default_output(fs, params, tailTotal);
    return;
end

% 2. 累积能量触发点
cumE = cumsum(envSm);
triggerIdx = find(cumE >= totalE * cumEnergyFrac, 1, 'first');
if isempty(triggerIdx), triggerIdx = 1; end

% 3. ✅ 修复：正确的物理延迟范围计算
% 基于全局信息估算合理延迟范围
% 方法1：使用互相关延迟，但限制范围
if delayCorr > 0
    % 限制delayCorr在合理范围内
    validDelayCorr = min(max(delayCorr, minPhysDelay), maxPhysDelay);
    expectedStart = max(minPhysDelay, validDelayCorr - delaySearchRadius);
    expectedEnd = min(searchEnd, validDelayCorr + delaySearchRadius);
else
    % 方法2：如果没有有效的互相关延迟，使用经验范围
    expectedStart = minPhysDelay;
    expectedEnd = min(searchEnd, maxPhysDelay);
end

% 确保范围有效且起始点小于结束点
if expectedStart >= expectedEnd
    if debugMode
        fprintf('  [DEB-range] 范围无效 [%d,%d]，使用默认范围\n', ...
            expectedStart, expectedEnd);
    end
    expectedStart = minPhysDelay;
    expectedEnd = min(searchEnd, maxPhysDelay);
end

% 如果范围仍然无效，使用触发点到搜索结束
if expectedStart >= expectedEnd
    expectedStart = triggerIdx;
    expectedEnd = searchEnd;
end

expectedDelayRange = [expectedStart, expectedEnd];

% 4. 在期望延迟范围内搜索峰值
searchRange = expectedStart:expectedEnd;
if isempty(searchRange)
    pkLocal = find(envSm(triggerIdx:searchEnd) > th, 1, 'first');
    if isempty(pkLocal)
        [~,pkLocal] = max(envSm(triggerIdx:searchEnd));
        pkLocal = pkLocal + triggerIdx - 1;
    else
        pkLocal = pkLocal + triggerIdx - 1;
    end
else
    % 在期望范围内搜索
    envSm_range = envSm(searchRange);
    th_range = th;
    
    pkInRange = find(envSm_range > th_range, 1, 'first');
    if isempty(pkInRange)
        % 在范围内没有超过阈值的点，取范围内的最大值
        [~, maxIdx] = max(envSm_range);
        pkLocal = expectedStart + maxIdx - 1;
    else
        pkLocal = expectedStart + pkInRange - 1;
    end
end

% 5. 二次验证：如果峰值太小，尝试在整个范围内搜索
if envSm(pkLocal) < th * 0.5
    [~, maxIdx] = max(envSm(triggerIdx:searchEnd));
    pkLocal = maxIdx + triggerIdx - 1;
end

%% 峰值细化
if getP(params, 'peakRefineEnable', true)
    peakRefineRadius = getP(params, 'peakRefineRadius', 100);
    if pkLocal > 1 && pkLocal < length(h_full)
        refineStart = max(1, pkLocal - peakRefineRadius);
        refineEnd = min(length(h_full), pkLocal + peakRefineRadius);
        refineSegment = h_full(refineStart:refineEnd);
        
        % 使用插值法精确峰值定位
        [maxVal, maxIdx] = max(abs(refineSegment));
        if maxIdx > 1 && maxIdx < length(refineSegment)
            % 三点抛物线插值
            y1 = abs(refineSegment(maxIdx-1));
            y2 = abs(refineSegment(maxIdx));
            y3 = abs(refineSegment(maxIdx+1));
            
            if y1 > 0 && y2 > 0 && y3 > 0
                delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
                refinedIdx = maxIdx + delta;
                pkLocal = refineStart - 1 + refinedIdx;
            end
        end
    end
end

%% ✅ 修复：改进的peakReliability判定
% 1. 计算峰值窗口
peakWinRadius = min(50, floor(length(h_full)/20));
peakWinStart = max(1, floor(pkLocal) - peakWinRadius);
peakWinEnd = min(length(h_full), ceil(pkLocal) + peakWinRadius);
peakWin = h_full(peakWinStart:peakWinEnd);

% 2. 计算两种比例
peakAbsFrac = sum(abs(peakWin)) / (sum(abs(h_full)) + 1e-12);
peakEnergyFrac = sum(peakWin.^2) / (sum(h_full.^2) + 1e-12);

% 3. ✅ 修复：改进的侧瓣抑制计算
% 先找到真正的峰值（可能不是pkLocal）
[~, truePeakIdx] = max(abs(h_full));
if abs(truePeakIdx - pkLocal) > 50
    % 如果检测到的峰值与最大峰值相差较大，使用最大峰值
    if debugMode
        fprintf('  [DEB-peak] 检测峰=%0.1f, 真实峰=%d, 使用真实峰\n', ...
            pkLocal, truePeakIdx);
    end
    pkLocal = truePeakIdx;
    % 重新计算峰值窗口
    peakWinStart = max(1, floor(pkLocal) - peakWinRadius);
    peakWinEnd = min(length(h_full), ceil(pkLocal) + peakWinRadius);
    peakWin = h_full(peakWinStart:peakWinEnd);
    peakAbsFrac = sum(abs(peakWin)) / (sum(abs(h_full)) + 1e-12);
    peakEnergyFrac = sum(peakWin.^2) / (sum(h_full.^2) + 1e-12);
end

% 侧瓣区域定义
side_lobe_radius = min(200, floor(length(h_full)/4));
side_start = max(1, floor(pkLocal) - side_lobe_radius);
side_end = min(length(h_full), ceil(pkLocal) + side_lobe_radius);

% 主峰区域（排除区域）
main_peak_radius = min(30, side_lobe_radius/5);
exclude_start = max(1, floor(pkLocal) - main_peak_radius);
exclude_end = min(length(h_full), ceil(pkLocal) + main_peak_radius);

% 构建侧瓣区域掩码
side_mask = false(1, side_end - side_start + 1);
if exclude_start >= side_start && exclude_end <= side_end
    exclude_local_start = exclude_start - side_start + 1;
    exclude_local_end = exclude_end - side_start + 1;
    side_mask(exclude_local_start:exclude_local_end) = true;
end

% 提取侧瓣
side_lobe_samples = h_full(side_start:side_end);
side_lobe_samples(side_mask) = [];

if ~isempty(side_lobe_samples)
    main_peak_val = max(abs(h_full(exclude_start:exclude_end)));
    max_side_lobe = max(abs(side_lobe_samples));
    if max_side_lobe > 0
        side_lobe_suppression_db = 20*log10(main_peak_val/(max_side_lobe + 1e-12));
    else
        side_lobe_suppression_db = Inf;
    end
else
    side_lobe_suppression_db = Inf;
end

% 4. 预回声检查
pre_peak_region = h_full(1:max(1, floor(pkLocal)-10));
if ~isempty(pre_peak_region)
    pre_peak_energy = sum(pre_peak_region.^2);
    total_energy = sum(h_full.^2);
    pre_echo_ratio = pre_peak_energy / (total_energy + 1e-12);
else
    pre_echo_ratio = 0;
end

% 5. 综合判定条件
abs_ok = (peakAbsFrac >= minPeakFrac);
energy_ok = (peakEnergyFrac >= minPeakFrac * 0.3);

% ✅ 修复：基于测量质量的动态判定
% 如果侧瓣抑制为负，说明检测到错误的峰值
if side_lobe_suppression_db < 0
    % 尝试重新检测峰值
    [~, newPeak] = max(abs(h_full));
    if abs(newPeak - pkLocal) > 20
        pkLocal = newPeak;
        % 重新计算相关参数
        peakWinStart = max(1, floor(pkLocal) - peakWinRadius);
        peakWinEnd = min(length(h_full), ceil(pkLocal) + peakWinRadius);
        peakWin = h_full(peakWinStart:peakWinEnd);
        peakAbsFrac = sum(abs(peakWin)) / (sum(abs(h_full)) + 1e-12);
        peakEnergyFrac = sum(peakWin.^2) / (sum(h_full.^2) + 1e-12);
        
        % 重新计算侧瓣抑制
        exclude_start = max(1, floor(pkLocal) - main_peak_radius);
        exclude_end = min(length(h_full), ceil(pkLocal) + main_peak_radius);
        
        side_mask = false(1, side_end - side_start + 1);
        if exclude_start >= side_start && exclude_end <= side_end
            exclude_local_start = exclude_start - side_start + 1;
            exclude_local_end = exclude_end - side_start + 1;
            side_mask(exclude_local_start:exclude_local_end) = true;
        end
        
        side_lobe_samples = h_full(side_start:side_end);
        side_lobe_samples(side_mask) = [];
        
        if ~isempty(side_lobe_samples)
            main_peak_val = max(abs(h_full(exclude_start:exclude_end)));
            max_side_lobe = max(abs(side_lobe_samples));
            if max_side_lobe > 0
                side_lobe_suppression_db = 20*log10(main_peak_val/(max_side_lobe + 1e-12));
            else
                side_lobe_suppression_db = Inf;
            end
        else
            side_lobe_suppression_db = Inf;
        end
    end
end

% 基于预回声程度的动态侧瓣抑制要求
if pre_echo_ratio > 0.15  % 预回声严重
    side_lobe_ok = (side_lobe_suppression_db > 8);
elseif pre_echo_ratio > 0.05  % 中等预回声
    side_lobe_ok = (side_lobe_suppression_db > 3);
else  % 预回声较小
    side_lobe_ok = (side_lobe_suppression_db > 1);
end

% 最终可靠性判定
peakReliability = (abs_ok || energy_ok) && side_lobe_ok;

% 特殊情况：如果能量比例特别高，放宽要求
if peakEnergyFrac > 0.15
    peakReliability = true;
end

% 检查峰值位置是否在合理范围内
if pkLocal < minPhysDelay || pkLocal > maxPhysDelay
    if debugMode
        fprintf('  [WARN] 峰值位置%0.1f超出物理延迟范围[%d, %d]\n', ...
            pkLocal, minPhysDelay, maxPhysDelay);
    end
    if pkLocal < minPhysDelay
        peakReliability = false;
    end
end

% 调试信息
if debugMode
    fprintf('  [DEB-deconv] pk=%0.1f (range=[%d,%d]), absFrac=%.4f, energyFrac=%.4f, sideSupp=%.1f dB, preEcho=%.3f -> reliable=%d\n', ...
        pkLocal, expectedDelayRange(1), expectedDelayRange(2), ...
        peakAbsFrac, peakEnergyFrac, side_lobe_suppression_db, pre_echo_ratio, peakReliability);
end

%% 截窗保留前导（生成最终输出 IR）
winStartLocal = max(floor(pkLocal) - preDelayKeep, 1);
winStopLocal  = min(winStartLocal + tailTotal - 1, length(h_full));
h_out = h_full(winStartLocal:winStopLocal);

peakIdxFinal = floor(pkLocal) - winStartLocal + 1;
if peakIdxFinal < 1, peakIdxFinal = 1; end

% SNR 估计
bodyL = max(1, peakIdxFinal - snrBodyRadius);
bodyR = min(length(h_out), peakIdxFinal + snrBodyRadius);
bodySlice = h_out(bodyL:bodyR);
tailSlice = h_out(max(length(h_out)-4*snrBodyRadius,1):end);
snrEst = 20*log10((rms(bodySlice)+1e-12)/(rms(tailSlice)+1e-12));

preEnergy = sum(abs(h_out(1:peakIdxFinal-1)));
postEnergy = sum(abs(h_out(peakIdxFinal:end)));
preEnergyFrac = preEnergy / (preEnergy + postEnergy + 1e-12);

% IR质量指标
ir_quality = struct();
ir_quality.preEchoRatio = pre_echo_ratio;
ir_quality.peakSharpness = max(abs(h_out)) / (rms(h_out) + 1e-12);
ir_quality.energyConcentration = peakEnergyFrac;
ir_quality.sideLobeSuppression = side_lobe_suppression_db;

out = struct();
out.h = h_out;
out.delayCorr = delayCorr;
out.startIdxGlobal = startIdx;
out.peakIdx = peakIdxFinal;
out.peakIdxZeroBased = peakIdxFinal - 1;
out.peakReliability = peakReliability;
out.peakAbsFrac = peakAbsFrac;
out.peakEnergyFrac = peakEnergyFrac;
out.sideLobeSuppression = side_lobe_suppression_db;
out.preEnergyFrac = preEnergyFrac;
out.snrEst = snrEst;
out.noiseStd = noiseStd;
out.noiseMAD = noiseMAD;
out.thresholdUsed = th;
out.triggerIdx = triggerIdx;
out.pkLocalGlobal = pkLocal;
out.warnEarly = warnEarly;
out.sampleRate = fs;
out.paramsUsed = params;
out.irQuality = ir_quality;
end

function val = getP(p, name, defaultVal)
    if isstruct(p) && isfield(p,name) && ~isempty(p.(name))
        val = p.(name);
    else
        val = defaultVal;
    end
end

function out = default_output(fs, params, tailTotal)
    out.h = zeros(tailTotal, 1);
    out.delayCorr = 0;
    out.startIdxGlobal = 1;
    out.peakIdx = 1;
    out.peakIdxZeroBased = 0;
    out.peakReliability = false;
    out.peakAbsFrac = 0;
    out.peakEnergyFrac = 0;
    out.sideLobeSuppression = 0;
    out.preEnergyFrac = 0;
    out.snrEst = -Inf;
    out.noiseStd = 0;
    out.noiseMAD = 0;
    out.thresholdUsed = 0;
    out.triggerIdx = 1;
    out.pkLocalGlobal = 1;
    out.warnEarly = false;
    out.sampleRate = fs;
    out.paramsUsed = params;
    out.irQuality = struct('preEchoRatio', 0, 'peakSharpness', 0, 'energyConcentration', 0, 'sideLobeSuppression', 0);
end