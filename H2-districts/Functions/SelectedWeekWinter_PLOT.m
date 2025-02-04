function [] = SelectedWeekWinter_PLOT(linew,font,Time,start,finish,P_PV_opt,P_FC_opt,P_imp_opt,P_exp_opt,P_FC_nom,SOC_opt,P_load,FC_On_opt)

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime([start_date:(1/24):finish_date], 'ConvertFrom', 'datenum');

colors = [1 44 86 129 172 214];
cmap = crameri('batlow');

figure
area(tt,-P_PV_opt(start:finish), 'FaceColor', cmap(colors(2), :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
h1=plot(tt,-P_PV_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(2), :));
hold on
area(tt, -P_FC_opt(start:finish), 'FaceColor', cmap(colors(3), :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
h2=plot(tt,-P_FC_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(3), :));
hold on
area(tt, -P_imp_opt(start:finish), 'FaceColor', cmap(colors(4), :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
h3=plot(tt,-P_imp_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(4), :))
hold on 
area(tt,P_exp_opt(start:finish), 'FaceColor', cmap(colors(5), :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on 
h4=plot(tt,P_exp_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(5), :));
hold on
area(tt, P_load(start:finish), 'FaceColor', cmap(colors(6), :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on
h5=plot(tt,P_load(start:finish),'LineWidth',linew,'Color',cmap(colors(6), :));
hold on
grid on

xlim([min(tt) max(tt)])
% ylim([-1.2*max([max(P_imp_opt) max(P_e_opt)]) 1.2*max(P_PV_opt) ])
ylabel ('Power [kW]','fontweight','bold');
xlabel('Time','fontweight','bold');
legend([h1, h2, h3, h4, h5],'P_{PV}','P_{PEM}','P_{imp}','P_{exp}','P_{load}','location','eastoutside')


set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', font-4); % Set axis font size and font name
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); % Set line width
title ('Power trend over a Selected Week in Winter')
% figure
% plot(tt,SOC_opt(start:finish),'-.','LineWidth',linew);
% hold on 
% grid on
% 
% ylabel ('SOC [-]','fontweight','bold');
% xlim([min(tt) max(tt)])
% ylim([0 1.2])
% 
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [0, 0, 7, 4.5]); % Width=3.5in, Height=2.5in
% set(gca, 'FontSize', font); % Set axis font size and font name
% set(findall(gcf, 'Type', 'line'), 'LineWidth', 1); % Set line width
figure
plot(tt, FC_On_opt(start:finish),tt, P_FC_opt(start:finish));
legend('On/Off status', 'Fuel Cell Power');
xlabel('Time (hours)');
ylabel('Power (kW)');
title ('WINTER WEEK')
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', font-4); % Set axis font size and font name
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); % Set line width
xlim([min(tt) max(tt)])
grid on
end