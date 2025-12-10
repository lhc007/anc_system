function ok = check_anc_dependencies(cfg)
% 检查 ANC 仿真依赖：
% - 输入噪声文件存在且可读
% - 反馈路径文件包含 F 且维度正确
% - 次级路径文件包含 secondary 结构体
% - 关键函数在路径上（提示）
fprintf('🔍 开始 ANC 依赖项检查...\n');
ok = true;

% 输入噪声
if ~isfile(cfg.inputAudioFile)
    fprintf('❌ 输入噪声文件不存在: %s\n', cfg.inputAudioFile); ok=false;
else
    try
        info = audioinfo(cfg.inputAudioFile);
        fprintf('✅ 输入音频：采样率 %d Hz，%d 通道，%d 样本\n', info.SampleRate, info.NumChannels, info.TotalSamples);
    catch ME
        fprintf('❌ 输入音频无法读取: %s\n', ME.message); ok=false;
    end
end

% 反馈路径
if ~isfile(cfg.feedbackPathFile)
    fprintf('❌ 反馈路径文件不存在: %s\n', cfg.feedbackPathFile); ok=false;
else
    S = whos('-file', cfg.feedbackPathFile);
    if ~any(strcmp({S.name}, 'F'))
        fprintf('❌ 反馈路径文件中缺少变量 F\n'); ok=false;
    else
        load(cfg.feedbackPathFile,'F');
        if ndims(F)~=3 || size(F,2)~=numel(cfg.micChannels.reference) || size(F,3)~=cfg.numSpeakers
            fprintf('❌ 反馈路径 F 尺寸错误，应为 [Lfb x %d x %d]\n', numel(cfg.micChannels.reference), cfg.numSpeakers); ok=false;
        else
            fprintf('✅ 反馈路径 F：尺寸 [%d x %d x %d]\n', size(F));
        end
    end
end

% 次级路径
if ~isfile(cfg.secondaryPathFile)
    fprintf('❌ 次级路径文件不存在: %s\n', cfg.secondaryPathFile); ok=false;
else
    S = whos('-file', cfg.secondaryPathFile);
    if ~any(strcmp({S.name}, 'secondary'))
        fprintf('❌ 次级路径文件中缺少 secondary 结构体\n'); ok=false;
    end
end

% 函数提示
needed = {'fxlms_recursive','anc_plot_results'};
for i=1:numel(needed)
    if exist(needed{i},'file')~=2
        fprintf('⚠️ 警告：函数 %s 不在路径上（相关功能不可用）\n', needed{i});
    end
end

if ok
    fprintf('🎉 所有依赖项检查通过！可以安全运行仿真。\n');
else
    fprintf('🛑 依赖项不完整，请修复后再运行。\n');
end
end