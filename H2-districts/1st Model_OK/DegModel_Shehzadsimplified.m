clear all
close all
clc
%% Input data

% Define the main path manually
    path_main = 'C:\Users\nanf\Desktop\H2-districts';
% path for input files
    path_input = fullfile(path_main, 'Input');
% Reading an Excel file from the Input directory
    excelFilePath = fullfile(path_input, 'DegradationModel');
    Input         = readtable(excelFilePath);
% Input data from paper
d_f            = 0.1; %0.02;          
P_FC_max       = 11.615;%12.1633; % [kW]
eff_FC         = 0.5;     % Nominal efficency
% Input vector from optimization problem
P_FC_opt       = Input.P_FC;                % optimal fuel cell power [kW]
FC_On_opt      = Input.FC_On;               % integer variable for fuel cell operation

% Time vector of simulation of 1 year
nHours         = 8760; 
Time           = (1: nHours)';
nH_lifetime  =  40000; %[h]
HHV         = 39.39 * 3.6 * 10^3;           % Hydrogen higher heating value [kJ/kg] 
HHV1        = 3.544;                        % [kWh/Nm3]
rho_H2      = 0.0898765;                    % H2 density @ T=0°C,1 bar [kg/m3]
% Initialization of degradation rate
fuel_cons        = zeros(nHours, 1);
fuel_cons(1)     = eff_FC*HHV1;  % Initial fuel consumption value [kWh/Nm^3] 

% Degradation calculation loop
for i = 2:nHours
    fuel_cons(i) = (1 - ((d_f /nH_lifetime) * FC_On_opt(i))) * fuel_cons(i-1); % [kWh/Nm3]
end
% Mass flow rate considering degradation
flow_H2     = ((P_FC_opt.*FC_On_opt)./ fuel_cons)* rho_H2;         % [kg/h]
flow_H2_nom = ((P_FC_opt.*FC_On_opt)./eff_FC/HHV1)* rho_H2;        % [kg/h]
%% No sense 
tot_m_flow_H2            = sum(flow_H2);     
tot_m_flow_H2_nom        = sum(flow_H2_nom);
increase_mflowH2         =(tot_m_flow_H2-tot_m_flow_H2_nom)*100/tot_m_flow_H2_nom; % [%]
% eff_FC_deg               = P_FC_opt.*FC_On_opt*3600./flow_H2./HHV;
eff_FC_deg1              = fuel_cons./ HHV1;
annualEffReduction       = abs(eff_FC_deg1(end)-eff_FC_deg1(1))*100/eff_FC_deg1(1); % percentage annual efficency reduction [%]
New_Lifetime       = 10*1/annualEffReduction;

% Plot of degradation rate over time
figure;
plot(Time, eff_FC_deg1, 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Efficnency Degradation Rate');
title('Degradation efficency rate with respect Operational Hours');
