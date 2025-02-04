close all
clear all
clc

%% INPUT DATA
font   = 18;
linew  = 1.5; 
nHours = 8760; 
Time =(1:nHours)';
% operating parameters
k_r  = 1; %=1.72  % Correction coefficient of real-world driving cycles
n11  = 1412/8760;%1626/8760;        % time duration of heavy loading [cycles/h]
n12  = 157/8760;%401/8760;          % number of start/stop cycles [cycles/h]
t11  = 25/8760;%29/8760;            % time duration of low-power operation [h/h]
t12  = 4631/8760;%4359/8760;        % time duration of high-power operation [h/h]

% Coefficient in degradation rate from manufacturer's datasheet of a PEMFC (paper)
k1 = 0.0000593;     % [%/cycle] Degradation coefficient for load changing
k2 = 0.00196;       % [%/cycle] Degradation coefficient for start/stop cycling
k3 = 0.00126;       % [%/h]     Degradation coefficient for low-power operation (idle time)
k4 = 0.00147;       % [%/h]     Degradation coefficient for high-power operation

%% Fuel Cell Degradation Rate

D_fc1  = (k1 * n11 + k2 * n12 + k3 * t11 + k4 *t12); % Degradation rate [%/h]
D_lowpow1    = k3*t11;  %[%/h]
D_highpow1   = k4*t12;  %[%/h]
D_startstop1 = k2*n12;  %[%/h]
D_loading1   = k1*n11;  %[%/h]
% Lifetime predicted
Newlifetime1 = 10/D_fc1; %[h]

%% plot
degcoeffs= [D_lowpow1, D_highpow1,D_loading1,D_startstop1];
labels = {'Low Power Operation', 'High Power Operation', 'Loading Operation', 'Start/Stop cycling'};
cmap = crameri('batlow', length(degcoeffs));
figure;
b = bar(degcoeffs);
b.FaceColor = 'flat';
b.FaceAlpha = 0.8;
for i = 1:length(degcoeffs)
    b.CData(i, :) = cmap(i, :);
end

set(gca, 'xticklabel', labels);
set(gca, 'FontSize', font-4);
ylabel('Degradation rate (%)','FontWeight','bold');
ylim([0 3.5e-3])

