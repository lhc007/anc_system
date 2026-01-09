function metrics = assess_ir_quality(ir, snrList, cohList, cfg)
% assess_ir_quality - 综合评估IR质量（兼容没有 nan* 函数的 MATLAB）
% ir: [L x M] 脉冲响应
% snrList: [R x M] SNR矩阵（可能含 NaN）
% cohList: [R x M] 相干性矩阵（可能含 NaN）
% cfg: 配置结构体

    numMics = size(ir, 2);
    L = size(ir, 1);  % IR长度

    % ----- nan-safe 全局统计（兼容不同 MATLAB 版本） -----
    try
        medianSNR = median(snrList(:), 'omitnan');
        meanCoherence = mean(cohList(:), 'omitnan');
    catch
        % 手动实现 omitnan
        v1 = snrList(:);
        v1 = v1(~isnan(v1));
        if isempty(v1), medianSNR = NaN; else medianSNR = median(v1); end

        v2 = cohList(:);
        v2 = v2(~isnan(v2));
        if isempty(v2), meanCoherence = NaN; else meanCoherence = mean(v2); end
    end

    % 检查每个通道的质量（nan-safe）
    channelOK = true(1, numMics);
    for ch = 1:numMics
        try
            chSNR = median(snrList(:, ch), 'omitnan');
            chCoh = mean(cohList(:, ch), 'omitnan');
        catch
            vch1 = snrList(:, ch);
            vch1 = vch1(~isnan(vch1));
            chSNR = isempty(vch1) && NaN || median(vch1);

            vch2 = cohList(:, ch);
            vch2 = vch2(~isnan(vch2));
            chCoh = isempty(vch2) && NaN || mean(vch2);
        end

        if isempty(chSNR) || isnan(chSNR) || isempty(chCoh) || isnan(chCoh)
            channelOK(ch) = false;
        elseif chSNR < cfg.snrThresholdDB || chCoh < cfg.coherenceThreshold
            channelOK(ch) = false;
        end
    end

    % IR自身质量检查：能量比
    ir_energy = sum(ir.^2, 1);
    if all(ir_energy > 0)
        energy_ratio = max(ir_energy) / min(ir_energy);
    else
        energy_ratio = Inf;
    end

    % 查找主峰位置 & 一致性
    if numMics > 1
        peak_positions = zeros(1, numMics);
        for ch = 1:numMics
            [~, peak_idx] = max(abs(ir(:, ch)));
            peak_positions(ch) = peak_idx;
        end
        peak_consistency = std(peak_positions) / (mean(peak_positions) + eps);
    else
        peak_consistency = 0;
    end

    % 对称性检查
    if L > 50
        half_len = floor(L/2);
        first_half_energy = sum(ir(1:half_len, :).^2, 1);
        second_half_energy = sum(ir(half_len+1:end, :).^2, 1);
        if all(first_half_energy > 0)
            symmetry_ratio = median(second_half_energy ./ first_half_energy);
        else
            symmetry_ratio = Inf;
        end
    else
        symmetry_ratio = 1;
    end

    % 输出
    metrics.medianSNR = medianSNR;
    metrics.meanCoherence = meanCoherence;
    metrics.channelOK = channelOK;
    metrics.allChannelsOK = all(channelOK);
    metrics.usable = metrics.allChannelsOK && ...
                     ~isnan(medianSNR) && medianSNR >= cfg.snrThresholdDB && ...
                     ~isnan(meanCoherence) && meanCoherence >= cfg.coherenceThreshold;
    metrics.energyRatio = energy_ratio;
    metrics.peakConsistency = peak_consistency;
    metrics.symmetryRatio = symmetry_ratio;
end