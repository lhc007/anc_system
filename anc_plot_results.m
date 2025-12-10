function anc_plot_results(logData)
% 绘制 ANC 结果：误差RMS、改善量、‖W‖、|W|max、步长、自适应开关
if isempty(logData.frame)
    fprintf('⚠️ logData 为空，无法绘图。\n'); return;
end

figure('Name','ANC 运行结果','Color','w');

subplot(3,1,1);
plot(logData.frame, logData.errRms,'b'); grid on;
ylabel('误差 RMS'); title('误差 RMS 与改善量');
yyaxis right; plot(logData.frame, logData.improv,'r'); ylabel('改善量 (dB)');

subplot(3,1,2);
plot(logData.frame, logData.Wnorm,'k'); hold on; grid on;
plot(logData.frame, logData.Wmax,'m--');
ylabel('权重范数/最大值');
legend('‖W‖_2','|W|max');

subplot(3,1,3);
plot(logData.frame, logData.mu,'g'); hold on; grid on;
stairs(logData.frame, double(logData.adaptEnable),'c');
ylabel('步长 / 自适应标志'); xlabel('帧');
legend('μ','AdaptEnable');
end