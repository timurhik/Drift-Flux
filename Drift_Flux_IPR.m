%
% Script file Drift_Flux IPR plots combined IPR and two intake curves for
% both Mukherjee and Brill and Drift-Flux correlation models

clear variables
close all

% ---------------------------------------------------------------------------------------------
% Input data:
% ---------------------------------------------------------------------------------------------
av         = 0; % use reservoir pressure at external boundary
beta       = 82e6 ; % Forcheimer coefficient, m^-1. Not relevant for fluid 1.
f_w_sc     = 0.0; % = q_w_sc/(q_w_sc+q_o_sc) = water cut, -. Not relevant for fluids 1 and 2
fluid      = 3; % fluid  = 1: single-phase oil flow
                % fluid  = 2: single-phase gas flow
                % fluid >= 3: multi-phase gas-oil-water flow
h          = 20.0; % reservoir height, m
k          = 4e-14; % effective permeability for fluids 1 and 2, or absolute permeability
                       % for fluid 3, m^2
n_pt       = 100; % number of points in plot, - 
oil = 1; % parameter to select black oil model or volatile oil table, -. Not relevant for
%          fluids 1 and 2
%          oil = 1: black oil; parameters computed with the aid of Standing correlations
%          oil = 2: black oil; parameters computed with the aid of Glaso correlations
%          oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
p_R        = 26.0e6 % reservoir pressure, Pa
q_o_sc_max = -3e-3; % maximum oil flow rate, m^3/s
r_e        = 800; % external radius, m
r_w        = 0.1778; % well bore radius, m
R_go       = 50; % producing GOR as observed at surface, m^3/m^3. Not relevant for fluid 2
rho_g_sc   = 0.95; % gas density at standard conditions, kg/m^3. 
rho_o_sc   = 850;  % oil density at st. conditions, kg/m^3. Not relevant for fluid 2
rho_w_sc   = 1000; % water density at st. conditions, kg/m^3. Not relevant for fluids 1 and 2
semi       = 0; % semi = 0: use steady-state conditions
%                 semi = 1: use semi-steady-state conditions
simp       = 1; % simp = 0: numerical solution using Runge Kutta integration
%                 simp = 1: simplified (semi-)analytical solution
S          = 0; % skin, -
T_R        = 120; % reservoir temperature, deg. C

% Relative permeability data (not relevant for fluids 1 and 2):
k_rg_0  = 0.7; % end-point relative permeability to gas, -
k_ro_0  = 0.9; % end-point relative permeability to oil, -
k_rw_0  = 0.5; % end-point relative permeability to water, -
n_g     = 3; % Corey exponent for gas, -
n_og    = 3; % Corey exponent for oil in gas-oil flow, -
n_ow    = 3; % Corey exponent for oil in oil-water flow, -
n_w     = 3; % Corey exponent for water, -
S_gc    = 0.05; % critical gas saturation, -
S_or    = 0.10; % residual oil saturation, -
S_wi    = 0.15; % immobile water saturation, -
% --------------------------------------------------------------------------------------------
% End of input data
% ---------------------------------------------------------------------------------------------

% Check input:
if (simp == 1 && beta ~= 0)
    warning('Forcheimer flow only available numerically.')
    warning('Parameter "simp" reset to 0.')
    simp = 0;
end
if (simp == 0 && semi == 1)
    warning('Semi-steady state solution only available analytically.')
    warning('Parameter "semi" reset to 0.')
    semi = 0;
end
if (simp == 0 && av == 1)
    warning('Average reservoir pressure solution only available analytically.')
    warning('Parameter "av" reset to 0.')
    av = 0;
end

% Create data vectors:
rel = [k_rg_0,k_ro_0,k_rw_0,n_g,n_og,n_ow,n_w,S_gc,S_or,S_wi];
rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc];

% Compute and print bubble point pressure (for info only):
if fluid ~= 2
    switch oil
        case 1
            p_b = pres_bub_Standing(R_go,rho_g_sc,rho_o_sc,T_R)
        case 2
            p_b = pres_bub_Glaso(R_go,rho_g_sc,rho_o_sc,T_R)
        case 3
            p_b = pres_bub_volatile_oil(T_R)
    end
end
    
% Compute and plot IPR:
switch fluid
       otherwise % multi-phase       
        delta_q_o_sc = q_o_sc_max/n_pt; % oil rate increment, m^3/s
        p_wf = p_R;
        results = zeros(n_pt,2);
        for i=1:n_pt
            q_o_sc = i * delta_q_o_sc % oil rate, m^3/s
            q_g_sc = R_go * q_o_sc; % gas rate, m^3/s
            q_w_sc = (f_w_sc/(1-f_w_sc)) * q_o_sc; % water rate, m^3/s
            q_sc = [q_g_sc,q_o_sc,q_w_sc];
            if simp == 0 % numerical solution
                p_old = p_wf;
                p_wf = res(beta,fluid,h,k,oil,p_R,q_sc,r_e,r_w,rel,rho_sc,T_R);
                Delta_p = p_old - p_wf;
                results(i,1:2)= [-q_o_sc p_wf];
                if p_wf < 1e5 + Delta_p % to avoid pressures below atmospheric
                    break
                end
            else % Vogel solution
                p_old = p_wf;
                p_wf = res_Vogel(av,h,k,oil,p_R,q_sc,r_e,r_w,rel,rho_sc,semi,S,T_R);
                Delta_p = p_old - p_wf;
                results(i,1:2)= [-q_o_sc p_wf];
                if p_wf < 1e5 + Delta_p % to avoid pressures below atmospheric
                    break
                end
            end
        end
        nonzero_rows = results(:,1) ~= 0; % select non-zero rows
        results = results(nonzero_rows,:); % remove zero rows
        plot(results(:,1)*1e3,results(:,2)/1e6,'LineWidth',1)
        xlabel('Oil Flowrate,\it -q_{o,sc}\rm (10^{-3} m^3/s)')
        ylabel('FBHP,\it p_{wf}\rm (MPa)')
        grid on
end
hold on
% Plot an intake curve for Drift-flux model
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
q_o_sc_max = -3e-3; % maximum oil rate, m^3/s. Not relevant for fluid 2
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

% Compute and plot intake curve:
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

hold on
% Plot an intake curve for Mukherjee and Brill correlation
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
        plot(results(:,1)*1e3,results(:,2)/1e6,'LineWidth',2)
        xlabel('Oil Flow Rate,\it -q_{o,sc}\rm (10^{-3} m^3/s)')
        ylabel('FBHP,\it p_{wf}\rm (MPa)')
        legend('IPR curve','Drift-Flux intake curve','Mukherjee and Brill intake curve')
        grid on
end