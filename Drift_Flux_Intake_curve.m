% The script file Drift_Flux_Intake_curve.m plots two intake curves for
% Mukherjee and Brill and Drift-Flux correlations

clear variables
close all

% ---------------------------------------------------------------------------------------------
% Input data:
% ---------------------------------------------------------------------------------------------
alpha = from_deg_to_rad(60); % wellbore inclination , rad
d = 0.062; % tubing inside diameter, m
e = 30e-6; % tubing roughness, m
f_w_sc = 0.0; % = q_w_sc/(q_w_sc+q_o_sc) = water cut, -. Not relevant for fluids 1 and 2
fluid = 6; % fluid = 1: single-phase oil flow
%            fluid = 2: single-phase gas flow
%            fluid = 3: multi-phase gas-oil-water flow, Hagedorn and Brown correlation
%            fluid = 4: multi-phase gas-oil-water flow, Mukherjee and Brill correlation
%            fluid = 5: multi-phase gas-oil-water flow, Beggs and Brill correlation
%            fluid = 6: multi-phase gas-oil-water flow, Shi et al. drift flux model (exercise)
n_pt = 50; % number of points in plot, - 
oil = 1; % parameter to select black oil model or volatile oil table, -. Not relevant for
%          fluids 1 and 2
%          oil = 1: black oil; parameters computed with the aid of Standing correlations
%          oil = 2: black oil; parameters computed with the aid of Glaso correlations
%          oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
%          Note: When using a volatile oil table, the input parameters R_sb and rho_sc should 
%          be consistent with those used to generate the tabulated values. For
%          vol_oil_table_01: rho_g_sc = 0.80 kg/m^3, rho_o_sc = 800 kg/m^3, R_sb = 450 m^3/m^3 
p_tf = 0.5e6; % FTHP, Pa
q_g_sc_max = -8; % maximum gas rate, m^3/s. Not relevant for fluids 1, 3, 4 and 5 
q_o_sc_max = -5e-3; % maximum oil rate, m^3/s. Not relevant for fluid 2
R_go = 50; % producing GOR as observed at surface, m^3/m^3. Not relevant for fluid 2 
rho_g_sc = 0.95; % gas density at standard conditions, kg/m^3 
rho_o_sc = 850;  % oil density at standard conditions, kg/m^3. Not relevant for fluid 2
rho_w_sc = 1050; % water density at standard conditions, kg/m^3. Not relev. for fluids 1 and 2
T_tf = 30; % tubing head temperature, deg. C
T_wf = 120; % bottomhole temperature, deg. C
z_tvd = 3000; % true-vertical depth, m
% ---------------------------------------------------------------------------------------------
% End of input data
% ---------------------------------------------------------------------------------------------

% Compute auxiliary variables:
s_tot = z_tvd/cos(alpha); % total along-hole well depth, m 

% Create data vector:
rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc];

% Compute and plot intake curve for Drift-Flux model:
switch fluid
          
    case {3,4,5,6} % multi-phase
        delta_q_o_sc = q_o_sc_max/n_pt; % oil rate increment, m^3/s
        results = zeros(n_pt,2);
        for i=1:n_pt
            q_o_sc = i * delta_q_o_sc; % oil rate, m^3/s
            q_g_sc = R_go * q_o_sc; % gas rate, m^3/s
            q_w_sc = (f_w_sc/(1-f_w_sc)) * q_o_sc; % water rate, m^3/s
            q_sc = [q_g_sc,q_o_sc,q_w_sc];
            p_wf = pipe(alpha,d,e,fluid,oil,p_tf,q_sc,rho_sc,0,s_tot,T_tf,T_wf);
            results(i,1:2)= [-q_o_sc p_wf];
        end
        plot(results(:,1)*1e3,results(:,2)/1e6,'LineWidth',1)
        xlabel('Oil Flow Rate,\it -q_{o,sc}\rm (10^{-3} m^3/s)')
        ylabel('FBHP,\it p_{wf}\rm (MPa)')
        grid on
end
% Adding the curve for MB model at the same plot
hold on
fluid = 4;
switch fluid
    case {3,4,5,6} % multi-phase
        delta_q_o_sc = q_o_sc_max/n_pt; % oil rate increment, m^3/s
        results = zeros(n_pt,2);
        for i=1:n_pt
            q_o_sc = i * delta_q_o_sc; % oil rate, m^3/s
            q_g_sc = R_go * q_o_sc; % gas rate, m^3/s
            q_w_sc = (f_w_sc/(1-f_w_sc)) * q_o_sc; % water rate, m^3/s
            q_sc = [q_g_sc,q_o_sc,q_w_sc];
            p_wf = pipe(alpha,d,e,fluid,oil,p_tf,q_sc,rho_sc,0,s_tot,T_tf,T_wf);
            results(i,1:2)= [-q_o_sc p_wf];
        end
        plot(results(:,1)*1e3,results(:,2)/1e6,'LineWidth',1)
        xlabel('Oil Flow Rate,\it -q_{o,sc}\rm (10^{-3} m^3/s)')
        ylabel('FBHP,\it p_{wf}\rm (MPa)')
        legend('Mukherjee and Brill correlation intake curve','Drift-Flux correlation intake curve') % Legend for all curves
        grid on
end

