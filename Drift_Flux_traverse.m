% The scipt file Drift_Flux_traverse plots pressure traverse for Drift-Flux
% model
clear variables
close all

% ---------------------------------------------------------------------------------------------
% Input data:
% ---------------------------------------------------------------------------------------------
alpha = from_deg_to_rad(60); % wellbore inclination , rad
d = 0.062; % tubing inside diameter, m 
e = 30e-6; % tubing roughness, m
fluid = 6; % fluid = 1: single-phase oil flow
%            fluid = 2: single-phase gas flow
%            fluid = 3: multi-phase gas-oil-water flow, Hagedorn and Brown correlation
%            fluid = 4: multi-phase gas-oil-water flow, Mukherjee and Brill correlation
%            fluid = 5: multi-phase gas-oil-water flow, Beggs and Brill correlation
%            fluid = 6: multi-phase gas-oil-water flow, Shi et al. drift flux model (exercise)
oil = 1; % parameter to select black oil model or volatile oil table, -. Not relevant for
%          fluids 1 and 2
%          oil = 1: black oil; parameters computed with the aid of Standing correlations
%          oil = 2: black oil; parameters computed with the aid of Glaso correlations
%          oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
%          Note: When using a volatile oil table, the input parameters R_sb and rho_sc should 
%          be consistent with those used to generate the tabulated values. For
%          vol_oil_table_01: rho_g_sc = 0.80 kg/m^3, rho_o_sc = 800 kg/m^3, R_sb = 450 m^3/m^3 
p_tf = 0.5e6 % FTHP, Pa 
R_go = 50; % Gas oil ratio, m3/m3
q_g_sc = -0.0001*R_go; % gas rate at standard conditions, m^3/s
q_o_sc = -0.0001; % oil rate at standard conditions, m^3/s. Should be 0 for fluid 2
q_w_sc = -0.2e-2; % water rate at standard conditions, m^3/s. Should be 0 for fluids 1 and 2
% Note: flow rates should have negative values for a production well!
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

% Create data vectors:
q_sc = [q_g_sc,q_o_sc,q_w_sc];
rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc];


% Compute and plot traverse top-down:
[p_wf,s,p] = pipe(alpha,d,e,fluid,oil,p_tf,q_sc,rho_sc,0,s_tot,T_tf,T_wf);
p_wf
figure
plot(p/1e6,s,'LineWidth',1)
axis ij
xlabel('Wellbore Pressure,\it p\rm (MPa)')
ylabel('AHD,\it s\rm (m)')
grid on
legend('\itp_{\rmtot}','\itp_{\rmgrav}','\itp_{\rmfric}','\itp_{\rmacc}','location','NorthEast')
