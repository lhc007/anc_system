clear; clc;

%% 1. 查询所有播放设备（OUTPUT）
fprintf('=== 播放设备（扬声器）列表 ===\n');

    info = audiodevinfo;
    if isfield(info, 'output')
        for i = 1:length(info.output)
            fprintf('设备ID: %d, 名称: %s\n', info.output(i).ID, info.output(i).Name);
        end
    end


fprintf('\n=== 录音设备（麦克风）列表 ===\n');
%% 2. 查询所有录音设备（INPUT）
    info = audiodevinfo;
    if isfield(info, 'input')
        for i = 1:length(info.input)
            fprintf('设备ID: %d, 名称: %s\n', info.input(i).ID, info.input(i).Name);
        end
    end
