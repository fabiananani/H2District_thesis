function [] = PowerFCSummer_cost_PLOT(linew,font,numCosts,Time,start,finish,P_FC_opt_matrix,S_FC_matrix,h2_cost)
% this function enables to plot the trend of the power and the optimal
% size of the Fuel Cell with respect to the different costs of H2 over a
% selected week 

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime([start_date:(1/24):finish_date], 'ConvertFrom', 'datenum');


% Define colors for each hydrogen price plot
cmap = crameri('batlow', numCosts); 
colors = cmap;  
lineStyles = {'-', '-', ':', '-.'};
%% plot with 2 y-axis
figure;
hold on;
grid on;
% Create primary y-axis for FC power
yyaxis left;
for i = 1:numCosts
        plot(tt, P_FC_opt_matrix(start:finish, i), 'Color', colors(i,:), ...
            'LineWidth',1,'LineStyle', lineStyles{mod(i-1, length(lineStyles)) + 1},'DisplayName', ['CostH2 = ' num2str(h2_cost(i)) ' CHF/kg',]);
end

ylabel('Power [kW]');
set(gca, 'YColor', 'k');  % Set left y-axis color to black

% xlabel('Date');
xlim([min(tt) max(tt)]);
% title('FC Power for Different Hydrogen Prices in Summer');

% Create secondary y-axis for S_FC
yyaxis right;
for i = 1:numCosts
    plot(tt, S_FC_matrix(start:finish, i), 'Color', colors(i,:), 'LineStyle', '--', ...
        'LineWidth',1,'DisplayName',['CostH2 = ' num2str(h2_cost(i)) ' CHF/kg']);
end

ylabel('Size [kW]');   
set(gca, 'YColor', [0.5, 0, 0.1]);

leg = legend('Location', 'eastoutside');
set(gca, 'FontSize', font-4); % Set axis font size and font name
set(gcf, 'Units', 'inches');
set(leg, 'FontSize',8);     
set(gcf, 'Position', [0, 0, 9 , 4.5]);
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); % Set line width
hold off;


end

