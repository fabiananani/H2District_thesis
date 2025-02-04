clear all
close all
clc

%% Set up Gurobi: Adding path to optimizer

addpath('C:\gurobi1103\win64\examples\matlab')
gurobi_setup
%% Input data

% Define the main path manually, e.g. 'C:\Users\My_models'
    path_main = 'C:\Users\nanf\Desktop\H2-districts';
% path for input files
    path_input = fullfile(path_main, 'Input');
% Define the path for functions
    path_functions = fullfile(path_main, 'Functions');

% Reading an Excel file from the Input directory
    excelFilePath1 = fullfile(path_input, 'Inputs_BoundaryLoads_2022.xlsx');
    excelFilePath2 = fullfile(path_input, 'nestel_demand.csv');
    excelFilePath3 = fullfile(path_input,'PriceImpExp_ele.xlsx');
    Input1         = readtable(excelFilePath1);
    Input2         = readtable(excelFilePath2);
    Input3         = readtable(excelFilePath3);
    
%% Pre-processing

    irradiance     = Input1.G;                     % Hourly solar irradiance    [W/m^2]
    T_amb          = Input1.T_amb;                 % Hourly temperature         [°C]
    nHours         = 8760;                        % Number of hours simulated (8760h)    
    idxHr2ToEnd    = (2:nHours)';                 % Hours until the end
    Time           = (1:nHours)';                 % Time vector
    days           = nHours/24;                   % Number of days simulated
    weeks          = days/7;                      % Number of weeks simulated
    clear t; 
    linew          = 1.5;             
    font           = 18;
    
%% INPUT PARAMETERS

 % general inputs

    bigM            = 60;      % Large number for big-M constraints
    MaxSimTime      = 900;     % Maximum time for MILP solver                 [s]      
    deltat          = 3600;    % time step di 1h                              [s]       

 % Electrical district demand
        
    P_load         = Input2.P_el_nest_h;       % NEST electrical demand      [kW]

 % Fixed technical parameters  

    HHV         = 39.39 * 3.6 * 10^3;           % Hydrogen higher heating value [kJ/kg]  
    
 % Flag for electricity cost (true/false)
   
    use_Empa_prices = true; % change with "false" to use DSO prices
    
    if use_Empa_prices

        c_gridimp   = Input3.PriceimpEmpa;                     % [CHF/kWh]
        c_gridexp   = Input3.PriceexpEmpa;                     % [CHF/kWh]
    else 
                c_gridimp   = Input3.Priceimpgreen;                  % [CHF/kWh]
                c_gridexp   = Input3.PriceexpPV;                     % [CHF/kWh]
    end 

    % Efficencies and max lîfe time from datasheet

    eff_b       = 0.95;                                  % constant efficiency of the battery (round trip efficency)
    eff_PV      = 0.17;%0.20;                            % constant efficiency of the PV system
    eff_FC      = 0.5;                                   % CONSTANT electrical efficency of FC  
    max_LT      = 40000;                                 % [h]

%range of H2 costs

    c_h2            = 4;    %[2 4 6 8];                         % cost of hydrogen [CHF/kg]
    % numCosts        = numel(c_h2);
    % P_FC_opt_matrix = zeros(nHours, numCosts);                  % Matrix to store FC power for each H2 cost
    % S_FC_opt_matrix = zeros(nHours, numCosts);
    

% Fuel Cell (FC)
    S_FC_min     = 0;                            % Minimum size FC             [kW]
    S_FC_max     = 60;                           % Maximum size FC             [kW]

% PV   
    P_PV_peak_ref  = 100;                 % Reference value peak power PV      [kW]
    Area_PV_ref = P_PV_peak_ref/eff_PV;   % Reference area PV                  [m2]
    Area_PV_min = 0;                      % Minimum area PV                    [m2]
    Area_PV_max = Area_PV_ref*2;          % Maximum area PV                    [m2]

% Battery

    C_b_max     = 96 * 3.6 * 10^3;   % Battery capacity                    [kJ]
    C_b_min     = C_b_max/3;         % MIN Battery capacity                [kJ]
    P_b_max     = C_b_max / 3600;    % Battery capacity (1C rate)          [kWh]

% unit prices and lifetime components

    d           = 0.05;%0.04;              % Discount rate
    ann         = d / (1 - (1 + d)^(-20));  % annuity factor calculated with plant lifetime        
   
    UP_PV       = 1300;%1240;         % Unit price PV                                    [CHF/kW_p]
    life_PV     = 25;%30;             % Lifetime PV                                      [years]
    maint_PV    = 0.01;%0.0158;       % Annual cost maintenance PV, frac total cost
    ann_PV      = d / (1 - (1 + d)^(-life_PV));
   
    UP_b        = 1000;%600      % Unit price battery [CHF/kWh]
    life_b      = 12;%10;        % Lifetime battery                                 [years]
    maint_b     = 0.02;          % frac capex 
    ann_b       = d / (1 - (1 + d)^(-life_b));
    
    UP_FC        = 950;%3000          % Unit price FC [CHF/kW]
    life_FC      = 7.5;%15            % Lifetime FC                                          [years]
    maint_FC     = 0.024;%0.02        % Annual cost maintenance FC, frac capex cost
    ann_FC       = d / (1 - (1 + d)^(-life_FC));
    
    % for i = 1:numCosts
    % Set hydrogen cost for this iteration
    current_c_h2 = c_h2;
%% Define the optimization problem and the optimization variables

    sizingprob = optimproblem;

%% DECISION VARIABLES

% SIZING decision variables

    % Battery capacity in [kJ]
    C_b            = optimvar('C_b','LowerBound',C_b_min,'UpperBound',C_b_max);  % C_b actual installed capacity,C_b_max= max capacity we can have
    % Fuel Cell size in [kW]
    S_FC           = optimvar('S_FC','LowerBound',S_FC_min,'UpperBound',S_FC_max);
    % PV area in [m2]
    Area_PV        = optimvar('Area_PV','LowerBound',Area_PV_min,'UpperBound',Area_PV_max);

% OPERATIONAL decision variables
    
    % fuel cell generated power in [kW]
    P_FC           = optimvar('P_FC',nHours,'LowerBound',0,'UpperBound',S_FC_max);
    % additional variables for FC operation implementation
    P_FC_On        = optimvar('P_FC_On',nHours,'LowerBound',0,'UpperBound',S_FC_max);
    delta_On          = optimvar('FC_On',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    % imported power from the grid in [kW]
    P_imp          = optimvar('P_imp',nHours,'LowerBound',0);
   
    % exported power to the grid in [kW]
    P_exp          = optimvar('P_exp',nHours,'LowerBound',0);
        
    % Battery energy content in [kWh]
    E_b            = optimvar('E_b',nHours,'LowerBound',0,'UpperBound',C_b_max);   
    % Battery charging power in [kW]
    P_b_ch         = optimvar('P_b_ch',nHours,'LowerBound',0,'UpperBound',P_b_max);
    % Battery discharging power in [kW]
    P_b_disch      = optimvar('P_b_disch',nHours,'LowerBound',0,'UpperBound',P_b_max);
    
    %% Derived variables

    % PV generated power
    P_PV      = irradiance.*eff_PV*Area_PV/1000;                           % [kW]
    P_PV_peak = 1000*eff_PV*Area_PV/1000;                                  % [kW]
    
    % C-rate of battery [kWh]
    P_b_lim     = C_b / 3600;                                              % new max capacity in [kWh], 1-C rate means that the battery charge/dsicharge in 1 hour

    % mass flow of H2
    m_flow_H2   = (P_FC/eff_FC/HHV)*3600;                                  % Consumed mass flow rate of H2 [kg/h]
    %% CONSTRAINTS

    % energy balance

    sizingprob.Constraints.EnBalance    = (P_FC + P_PV + P_b_disch + P_imp) == (P_b_ch + P_load + P_exp);
 
   
    % BATTERY

    % sizingprob.Constraints.NoSimultaneousChDisch  = discharging_on + charging_on <= ones(nHours,1);
    sizingprob.Constraints.PowerBatt_ch_0         = P_b_ch(1) == 0; 
    sizingprob.Constraints.PowerBatt_disch_0      = P_b_disch(1) == 0; 

    %sizingprob.Constraints.Ch_on1                = P_b_ch <= P_b_max * charging_on;
    sizingprob.Constraints.Ch_on2                 = P_b_ch <= P_b_lim;
    %sizingprob.Constraints.Disch_on1             = P_b_disch <= P_b_max * discharging_on;
    sizingprob.Constraints.Disch_on2              = P_b_disch <= P_b_lim;
    sizingprob.Constraints.E_b                    = E_b(idxHr2ToEnd) - E_b(idxHr2ToEnd-1) == P_b_ch(idxHr2ToEnd)*deltat - P_b_disch(idxHr2ToEnd)*deltat;
    sizingprob.Constraints.E_b_cont               = E_b(1) == E_b(end); 
    sizingprob.Constraints.E_b_max                = E_b <= C_b;
    sizingprob.Constraints.E_b_min                = E_b >= 0.2 * C_b;
    
    % FC    
    sizingprob.Constraints.MaxPowerFC             = P_FC <= P_FC_On;
    % sizingprob.Constraints.MinPowerFC             = P_FC >= 0.2 * P_FC_On;
    sizingprob.Constraints.PowerFC1               = P_FC_On <= S_FC_max * delta_On;
    sizingprob.Constraints.PowerFC2               = P_FC_On <= S_FC;
    sizingprob.Constraints.PowerFC3               = P_FC_On >= S_FC - S_FC_max.*(ones(nHours,1) - delta_On);
   
   
%% OBJECTIVE FUNCTION

cost_inst     = (P_PV_peak * UP_PV * ann_PV + S_FC * UP_FC * ann_FC + C_b/3600 * UP_b * ann_b)/1000; % [kEUR/y]
%cost_op      = sum(c_h2*m_flow_H2)/1000 + sum(P_imp.*c_gridimp)/1000 - sum(P_exp.*c_gridexp)/1000;  % [kEUR/y]
cost_imp      = sum(current_c_h2 * m_flow_H2) / 1000 + sum(P_imp .* c_gridimp) / 1000;               % [kEUR/y]
cost_exp      = sum(P_exp .* c_gridexp) / 1000;                                                      % [kEUR/y]
cost_maint    = (maint_PV * P_PV_peak * UP_PV + maint_FC * S_FC * UP_FC + maint_b * UP_b * C_b/3600)/1000;  % [kEUR/y]

cost = cost_inst + cost_imp - cost_exp + cost_maint;

% set objective
sizingprob.Objective = cost;

%% Solve optimization problem

intcon = [];
options = optimoptions('intlinprog','MaxTime',MaxSimTime);
[solution,fval,reasonSolverStopped] = solve(sizingprob,'Options',options);

% show problem
% show(sizingprob);

% % Store the results for this H2 costs
% P_FC_opt_matrix(:, i)   = solution.P_FC;  % Store FC power vector for this hydrogen cost
% S_FC_opt_matrix(:, i)   = solution.S_FC;
%     end
%% Post-processing and results overview

cost_total      = evaluate(cost, solution);                                      % optimization problem result
cost_inst       = evaluate(cost_inst, solution);                                 % kEUR/y
cost_imp        = evaluate(cost_imp, solution);                                  % kEUR/y
cost_maint      = evaluate(cost_maint, solution);                                % kEUR/y
cost_exp        = evaluate(cost_exp,solution);

Area_PV_opt     = solution.Area_PV;                                              % [m2]
P_PV_opt        = irradiance.*eff_PV.*solution.Area_PV./1000;                    % [kW]
P_FC_opt        = solution.P_FC;                                                 % [kW]
P_b_disch_opt   = solution.P_b_disch;                                            % [kW]
P_b_ch_opt      = solution.P_b_ch;                                               % [kW]
P_imp_opt       = solution.P_imp;                                                % [kW]
P_exp_opt       = solution.P_exp;                                                % [kW]
SOC_opt         = solution.E_b/solution.C_b;                                     % [-]
S_FC_opt        = solution.S_FC;                                                 % [kW]
C_b_opt         = solution.C_b;                                                  % [kJ]
E_b_opt         = solution.E_b;
E_ch            = sum(solution.P_b_ch);
E_FC            = sum(solution.P_FC);                                            %[kWh]

E_exp           = sum(solution.P_exp);
E_PV            = sum(P_PV_opt);
E_disch         = sum(solution.P_b_disch);
E_imp           = sum(solution.P_imp);
E_load          = sum(P_load); 
E_consumed      = E_ch + E_exp + E_load;
E_supplied      = E_PV + E_disch + E_imp + E_FC;

delta_On_opt    = evaluate(delta_On,solution); 
w_hours_opt     = sum(delta_On_opt);    

%% Start-Stop cycle per hour calculation
state_changes = diff(delta_On_opt);
% Identify start transitions (0 → 1) and stop transitions (1 → 0)
starts = state_changes == 1; % Logical array for start transitions
stops = state_changes == -1; % Logical array for stop transitions
% Total number of transitions (both starts and stops)
total_transitions = sum(starts) + sum(stops);
% Each start + stop pair constitutes one cycle
total_cycles = total_transitions / 2;

%% Time duration of low-power operation condition In Hours 
% Low power threshold (5% of maximum power)
low_power_threshold = 0.05 * S_FC_opt;
low_power_hours = (0 < P_FC_opt)&(P_FC_opt <= low_power_threshold);
% total duration in hours
low_power_duration = sum(low_power_hours); %[h]
%% Time duration of high-power operation condition in Hours
% High power threshold (90% of maximum power)
high_power_threshold = 0.9 * S_FC_opt;
high_power_hours = P_FC_opt >= high_power_threshold;
% total duration in hours
high_power_duration = sum(high_power_hours); %[h]
%% Time duration of large load variations
% Load variations threshold (10% of maximum power)
load_var_threshold = 0.10 * S_FC_opt;
load_var_rate      = diff(P_FC_opt)/1; % [kW]
load_cycles        = abs(load_var_rate) >= load_var_threshold;
num_load_cycles    = sum(load_cycles); % Count the number of cycles

%%  plots
costs = [cost_inst, cost_imp,-cost_exp,cost_maint,cost_total];
labels = {'Installation', 'Import', 'Export', 'Maintenance','Total'};
cmap = crameri('batlow', length(costs));
figure;
b = bar(costs);
b.FaceColor = 'flat';
b.FaceAlpha = 0.5;
for i = 1:length(costs)
    b.CData(i, :) = cmap(i, :);
end
grid on
set(gca, 'xticklabel', labels);
set(gca, 'FontSize', font-5);
% ylim([-cost_total cost_total+3])
ylabel('Costs (kEUR/y)');


% cmap = crameri('bamako', numCosts);
% colors = cmap;  
% figure('Position', [100, 100, 1700, 300]); 
% plot(Time, P_load,'Color',colors(4,:),'LineWidth',0.6)
% grid on
% xlabel('Time [h]','FontWeight', 'bold');
% xlim([Time(1) Time(end)]);
% ylabel('Pload [kW]','FontWeight', 'bold');
% 
% 
% figure('Position', [100, 100, 1700, 300]); 
% plot(Time, irradiance/1000,'Color',colors(3,:),'LineWidth',0.6)
% grid on
% xlabel('Time [h]','FontWeight', 'bold');
% xlim([Time(1) Time(end)]);
% ylabel('G [kW/m^2]','FontWeight', 'bold');

% %title ('Total Irradiance','FontWeight', 'bold')
% 
% figure('Position', [100, 100, 1700, 300]); 
% plot(Time(1:168),c_gridimp(1:168),'Color',colors(1,:),'LineWidth',linew)
% grid on
% xlabel('Time [h]','FontWeight', 'bold');
% xlim([Time(1) Time(168)]);
% ylabel('Cost [CHF/kWh]','FontWeight', 'bold');
% % title('Cost of imported electricity')
% 
% figure('Position', [100, 100, 1700, 300]); 
% plot(Time,c_gridexp,'Color',colors(2,:),'LineWidth',linew)
% grid on
% xlabel('Time [h]','FontWeight', 'bold');
% xlim([Time(1) Time(end)]);
% ylabel('Cost [CHF/kWh]','FontWeight', 'bold');
% % title('Cost of exported electricity','FontWeight', 'bold')

%% Figures - single week plots
 
% functions directory
addpath(path_functions);

% for a specific week between start and finish
startS=6*26*24;
finishS=startS+24*7+1;
startW=2*26*24;
finishW=startW+24*7+1;

SelectedWeekSummer_PLOT(linew,font,Time,startS,finishS,P_PV_opt,P_FC_opt,P_imp_opt,P_exp_opt,S_FC,SOC_opt,P_load,delta_On_opt)
movegui('center');
SelectedWeekWinter_PLOT(linew,font,Time,startW,finishW,P_PV_opt,P_FC_opt,P_imp_opt,P_exp_opt,S_FC,SOC_opt,P_load,delta_On_opt)
movegui('center');
% PowerFCSummer_cost_PLOT(linew,font,numCosts,Time,startS,finishS,P_FC_opt_matrix,S_FC_opt_matrix,c_h2)
% movegui('center');
% PowerFCWinter_cost_PLOT(linew,font,numCosts,Time,startW,finishW,P_FC_opt_matrix,S_FC_opt_matrix,c_h2)
% movegui('center');