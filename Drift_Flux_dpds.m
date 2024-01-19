function dpds = Drift_Flux_dpds(s,p,~,alpha,d,e,oil,q_sc,rho_sc,s_in,s_out,T_in,T_out)
% function dpds = Hag_Brown_dpds(s,p,~,alpha,d,e,oil,q_sc,rho_sc,s_in,s_out,T_in,T_out)   
%
% Computes the derivative dp/ds for a given pressure p and along-hole distance s, in an element
% of a flowline-wellbore system. The distance s is measured from the separator towards the
% reservoir. Therefore, flowrates are negative for production wells.
%
% Uses the Hagedorn and Brown correlation for multiphase flow in (near-)vertical wells; see
% references [1] and [2]. A reality check has been added to ensure that the computed liquid
% hold-up (for flow with slip) is never smaller than the in-situ liquid volume fraction (the
% 'no-slip hold-up'). 
% 
% The vector p contains the total pressure, and the pressures taking into account the
% individual effects of gravity, friction and acceleration losses respectively. Accordingly,
% the vector dpds contains the total pressure loss per unit length, as well as the individual
% gravity losses, friction losses and acceleration losses.  
%
% This function can be used to compute the pressure drop through numerical integration. It has
% the correct format to be used in conjunction with one of the standard numerical integration
% routines in MATLAB.   
%
% alpha = inclination wrt. vertical, rad; alternatively alpha can be a survey file (matrix)
%         with AHD values in the first column (in m) and inclination values in the second
%         column (in rad).   
% d = inside diameter, m
% dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc]
% dpds_acc  =  pressure gradient due to acceleration losses, Pa/m
% dpds_fric =  pressure gradient due to friction losses, Pa/m 
% dpds_grav = pressure gradient due to head losses, Pa/m 
% dpds_tot  = dpds_grav + dpds_fric + dpds_acc = total pressure gradient, Pa/m  
% e = roughness, m
% oil = parameter to select black oil model or volatile oil table, -  
%   oil = 1: black oil; parameters computed with the aid of Standing correlations
%   oil = 2: black oil; parameters computed with the aid of Glaso correlations
%   oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
% p = [p_tot,p_grav,p_fric,p_acc], Pa
% p_acc  = p_in + pressure increase (decrease for production wells) due to accel. losses, Pa
% p_fric = p_in + pressure increase (decrease for production wells) due to frict. losses, Pa
% p_grav = p_in + pressure increase (decrease for production wells) due to head loss, Pa
% p_tot  = p_in + pressure increase (decrease for production wells) due to gravity, friction
%          and acceleration losses, Pa 
% p_in = pressure at s_in, Pa 
% p_out = pressure at s_out, Pa
% q_sc = [q_g_sc,q_o_sc,q_w_sc], m^3/s
% q_g_sc = gas flow rate at standard conditions, m^3/s
% q_o_sc = oil flow rate at standard conditions, m^3/s
% q_w_sc = water flow rate at standard conditions, m^3/s
% rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc], kg/m^3 
% rho_g_sc = gas density at standard conditions, kg/m^3
% rho_o_sc = oil density at standard conditions, kg/m^3
% rho_w_sc = water density at standard conditions, kg/m^3
% s = along-hole distance, measured from the separator to the reservoir, m   
% s_in = starting point for the integration
% s_out = end point for the integration
% T_in = temperature at s_in, deg. C 
% T_out = temperature at s_out, deg. C 
%
% JDJ, 18-04-11, last revised 10-05-13
%
% References:
% [1] Hagedorn, A.R. and Brown, K.E., 1965: Experimental study of pressure gradients occurring
%     during continuous two-phase flow in small-diameter vertical conduits. Journal of
%     Petroleum Technology 17 (4) 475-484.   
% [2] Brill, J.P. and Mukherjee, H., 1999: Multiphase flow in wells, SPE Monograph Series,
%     vol.17, SPE, Richardson. 

% Compute internal variables:
epsilon = e/d; % dimensionless pipe roughness, -
g = 9.81; % acceleration of gravity, m/s^2

% Check sign of pressure:
p_tot = p(1); % first element of vector p is the total wellbore pressure, Pa 
if min(p_tot) < 1e5
    warning('Pressure below atmospheric.')
    dpds_tot  = 0; 
    dpds_grav = 0; 
    dpds_fric = 0; 
    dpds_acc  = 0; 
    dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc];
    return
end

% Determine inclination in case of survey file input:
if length(alpha) > 1
    n_sur = length(alpha(:,1)); % number of survey points
    if s < alpha(1,1)
        help = alpha(1,2);
    else if s > alpha(n_sur,1)
        help = alpha(n_sur,2);
        else
            help = interp1(alpha(:,1),alpha(:,2),s);
        end
    end
    clear alpha;
    alpha = help; % replace survey file by single inclination value, rad
end

% Compute local gas and liquid properties:
[mu_g,mu_l,q_g,q_l,rho_g,rho_l,sigma_gl,v_sg,v_sl] = local_gas_liq_props(d,oil,p_tot,...
    q_sc,rho_sc,s,s_in,s_out,T_in,T_out);


% Compute hold-ups (slip) and in-situ volume fractions (no-slip):
[H_g,H_l] = Drift_Flux_Holdup(alpha,d,q_g,q_l,rho_l,rho_g,sigma_gl)

% Check for free gas:
if abs(q_g) < 1.e-12 % no free gas - liquid flow only

    % Compute pressure gradient for liquid-only flow:
    v_l = v_sl; % local liquid velocity, m/s
    N_Re = rho_l*d*abs(v_l)/mu_l; % Reynolds number, -
    f = Moody_friction_factor(epsilon,N_Re); % friction factor, -
    dpds_grav = rho_l*g*cos(alpha); % gravity losses, Pa/m
    dpds_fric = -rho_l*f*v_l*abs(v_l)/(2*d); % friction losses, Pa/m
    dpds_acc = 0; % acceleration losses are neglegible, Pa/m
    dpds_tot = dpds_grav + dpds_fric + dpds_acc; % total pressure gradient, Pa/m
    dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc];
    
else % gas-liquid flow  

    % Compute 'slip' and 'no-slip' gas-liquid mixture properties:
    mu_ms = mu_l^H_l * mu_g^H_g; % 'no-slip' gas-liquid mixture viscosity (unusual definition),
    %                              Pa s 
    rho_ms = rho_l*H_l + rho_g*H_g; % 'slip' gas-liquid mixture density, kg/m^3
    v_ms = v_sg + v_sl; % local mixture velocity, m/s    
    
    % Compute pressure drop: 
    N_Re = rho_ms*abs(v_ms)*d/mu_ms; % gas-liquid mixture Reynolds number, -
    f = Moody_friction_factor(epsilon,N_Re); % gas-liquid mixture friction factor, -
    help21 = rho_ms*v_ms*v_sg/p_tot; % acceleration loss factor E_k, -
    if help21 >= 1 
        error('Choked flow. Reduce rate, increase diameter, or increase back pressure.')
    end
    help22 = rho_ms*g*cos(alpha);
    help23 = -f*rho_ms*v_ms*abs(v_ms)/(2*d); % alternative friction equation; different from
    %                                          original paper
    dpds_grav = help22; % gravity losses, Pa/m
    dpds_fric = help23; % friction losses, Pa/m
    dpds_acc = (help22+help23)*(help21/(1-help21)); % acceleration losses, Pa/m
    dpds_tot = dpds_grav + dpds_fric + dpds_acc; % total pressure gradient, Pa/m
    dpds = [dpds_tot;dpds_grav;dpds_fric;dpds_acc];
end

