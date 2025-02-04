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

% Input vector from optimization problem
P_FC_max       = 11.615;%12.21694;                  % [kW]
P_FC_opt       = Input.P_FC;                % optimal fuel cell power [kW]
FC_On_opt      = Input.FC_On;               % integer variable for fuel cell operation
% P_FC_years = repmat(P_FC_opt,5,1);
% FC_ON_years = repmat(FC_On_opt,5,1);
 
% Time vector of simulation of 1 year
font           = 18;
linew          = 1.5; 
nHours         = 8760; 
Time           = (1: nHours)';
eff_FC         = 0.5;                          % Nominal electrical efficency of FC 
HHV            = 39.39 * 3.6 * 10^3;           % Hydrogen higher heating value [kJ/kg]

% Lifetime in hours for PEMFC from paper of Shehzad et al. 
nH_lifetime  =  40000;                                  % [h]
% PEM fuel cell lifetime was defined by the time elapsed until 10% of the
% initial voltage/performance is lost fro STAYERS project
Delta_eff_FC = 0.1*0.5;                                 % [-]
% Fuel Cell degradation rate as a constant value in Nani-Model
Const = Delta_eff_FC / nH_lifetime;                     % [1/h]

%% Dynamic calculation for efficency to account for degradation based on state of operation ON/OFF

% Hourly values of degradated efficency starting from nominal value 
eff_deg_FC1     = zeros(nHours,1);
eff_deg_FC1(1)  = eff_FC;

% Cumulative summation of the operating hours (FC_On=1 if FC is ON, FC_On=0
% if FC is OFF) 
cumsum = zeros(nHours,1);
cumsum(1) = FC_On_opt(1);

% H2 consumed mass flow rate with respect a degradated efficency 

m_flow_H2 = (P_FC_opt.*FC_On_opt./eff_FC/HHV)*3600;  % "Nominal value": mass flow rate vector in no-degradation model 

for i=2:nHours
    cumsum(i)        =cumsum(i-1)+FC_On_opt(i);
    eff_deg_FC1(i)    = eff_FC-(Const*cumsum(i));  % original equation of the model for efficnecy degradation : eta_fc(t)= 0.5-CONST*cumsum(delta_ON), NOT IMPLEMENTED IN MILP
end 

 m_flow_H2new = (P_FC_opt./eff_deg_FC1./HHV)*3600; % massflowrate considering original equation for efficency degradation 
%% Post-processing
annualEffReduction1 = abs(eff_deg_FC1(end)-eff_deg_FC1(1))*100/eff_deg_FC1(1); % efficency reduction [%]
New_Lifetime1                = 10*1/annualEffReduction1;

tot_m_flow_H2_nom            = sum(m_flow_H2);     
tot_m_flow_H2_new1           = sum(m_flow_H2new);
increase_mflow1              =(tot_m_flow_H2_new1-tot_m_flow_H2_nom)*100/tot_m_flow_H2_nom; % [%]

%% plot
cmap = crameri('batlow',5);
colors = cmap;

%1
figure('Position', [100, 100, 800, 500]);
plot(Time,eff_deg_FC1,'LineWidth',3,'Color',colors(1,:))
ylabel('FC efficiency [-]','FontWeight','bold')
xlabel('Time [h]','FontWeight','bold')

%2
%period 
startS=6*26*24;
finishS=startS+24*7+1;
year=365*2023+126;
start_date=startS/24+year;
finish_date=finishS/24+year;
tt=datetime(start_date:(1/24):finish_date, 'ConvertFrom', 'datenum');

%figure
figure('Position', [100, 100, 1800, 300]);
plot(tt,m_flow_H2new(startS:finishS),'LineWidth',linew,'Color',colors(4,:))
hold on 
plot(tt,m_flow_H2(startS:finishS),'LineWidth',linew,'Color',colors(2,:))
hold on
legend('H2 Mass flow rate with degradation','H2 Mass flow rate without degradation')
ylabel('massflow rate [kg/h]','FontWeight','bold')
xlim([tt(1) tt(end)])
set(gca, 'FontSize', font-5);
leg = legend('Location', 'eastoutside');