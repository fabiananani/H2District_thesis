clear all
close all
clc
%% VERIFICATION OF PAPERS' MODEL 
% PAPER 1: "Lifetime prediction and the economic lifetime of Proton Exchange
% Membrane fuel cells" (2015)
% Fuel cell bus: PEMFC 
% P_rated = 10; % kW
% I_rated = 100;% A
% V_rated = P_rated*1e3/I_rated; % V
% n_cell = 100; 
% % Operating parameters per hour
% nn1 = 1.18; % cycle/h -start/stop cycle
% tt1 = 9/60; %h - average idling time per hour
% nn2 = 28; %cycle/h - average load change cycles per hour
% tt2 = 21/60; %h - average high power load operation time per hour
% % Voltage degradation rate 
% V1 = 13.79e-6; % V/cycle - START-STOP
% U1 = 8.662e-6; % V/h - IDLING
% V2 = 0.4185e-6; % V/cycle - LOAD CHANGE
% U2 = 10.00e-6; % V/h - HIGH POWER LOAD
% kk=1.72;
% % VOLTAGE DEGRADATION RATE 
% R1 = (nn1*V1+tt1*U1+nn2*V2+tt2*U2); % V/h
% T_f1  = 0.1/(1.72*R1);

%% 
% PAPER 2: "A quick evaluating method for automotive fuel cell lifetime" (2008)
% Coefficient in degradation rate from manufacturer's datasheet of a PEMFC (paper)
k1 = 0.0000593;     % [%/cycle] Degradation coefficient for load changing 
k2 = 0.00196;       % [%/cycle] Degradation coefficient for start/stop operation 
k3 = 0.00126;       % [%/h]     Degradation coefficient for low-power operation (idle time)
k4 = 0.00147;       % [%/h]     Degradation coefficient for high-power operation

% Operating parameters per hour
n21 = 56; %cycle/h - average load change cycles per hour
n22 = 0.99; % cycle/h -start/stop cycle
t21 = 13/60; %h/h - average idling time per hour
t22 = 14/60; %h/h - average high power load operation time per hour
% DEGRADATION RATE 
D_fc2 = (n21*k1+n22*k2+t21*k3+t22*k4); % [%/h]
Newlifetime2  = 10/D_fc2;
D_lowpow2    = k3*t21;  %[%/h]
D_highpow2   = k4*t22;  %[%/h]
D_startstop2 = k2*n22;  %[%/h]
D_loading2   = k1*n21;  %[%/h]
%% plot
degcoeffs= [D_lowpow2, D_highpow2,D_loading2,D_startstop2];
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
ylabel('Degradation rate (%)','FontWeight','bold');
ylim([0 3.5e-3])

%% 
% % PAPER 3: "Multi-objective optimization of hybrid PEMFC/Liion battery propulsion systems for small and
% %medium size ferries" (2021) 
% % DATA FROM PAPER 
% deltav_load  =0.0441; % [uV/kW]
% deltav_stup  = 23.91; % [uV/cycle]
% k1dv         = 0.0351; % [uV/A]
% k2dv         = 0.75; % [uV]
% k1i          = 4.3479; % [A/kW]
% k2i          = -39.729; % [A]
% % Simulation time 
% nHours = 8760; 
% Time =(1:nHours)';
% % Fuel cell maximum power (output of MILP)
% S_FC           = 12.21694;                  % [kW]
% % Input vector from optimization problem
% P_FC_opt       = Input.P_FC;                % optimal fuel cell power [kW]
% FC_On_opt      = Input.FC_On;               % integer variable for fuel cell operation
% % Nominal value of efficnency 
% eff_FC         = 0.5; 
% 
% %% COMPUTATION OF VOLTAGE DECREASE 
% dV_load=zeros(nHours,1);
% dV_load(1)= abs(P_FC_opt(1))*deltav_load;  %[uV]
% dV_stup=zeros(nHours,1);
% dV_stup(1)=FC_On_opt(1)*deltav_stup;
% Ifc=zeros(nHours,1);
% Ifc(1)  = k1i*P_FC_opt(1) + k2i*FC_On_opt(1); %[A]
% dV_Pfc=zeros(nHours,1);
% dV_Pfc(1)  = k1dv* Ifc(1)+k2dv; %[uV]
% dV_tot = zeros(nHours,1);
% dV_tot(1)=dV_load(1)+ dV_stup(1)+ dV_Pfc(1); %[uV]
% for i=2:nHours
% dV_load (i)= abs(P_FC_opt(i-1)-P_FC_opt(i))*deltav_load;  %[uV]
% dV_stup(i) = FC_On_opt(i)*deltav_stup; %[uV]
% Ifc(i)  = k1i*P_FC_opt(i) + k2i*FC_On_opt(i); %[A]
% dV_Pfc(i)  = k1dv* Ifc(i)+k2dv; %[uV]
% dV_tot(i)= dV_load(i)+ dV_stup(i)+ dV_Pfc(i); %[uV]
% end 
% DELTAV = sum(dV_tot)*1e-6; %[uV/year]
% 
