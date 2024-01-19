function [p_out,s,p] = pipe(alpha,d,e,fluid,oil,p_in,q_sc,rho_sc,s_in,s_out,T_in,T_out)
% function [p_out,s,p] = pipe(alpha,d,e,fluid,oil,p_in,q_sc,rho_sc,s_in,s_out,T_in,T_out) 
%
% Computes the pressure p_out at along-hole distance s_out in a deviated pipe element for a
% given pressure p_in at along-hole distance p_in, through numerical integration from s_in to
% s_out.  
%
% The function can be used in the following fluid modes:
%
% fluid = 1: single-phase oil flow
% fluid = 2: single-phase gas flow 
% fluid = 3: multi-phase gas-oil-water flow, Hagedorn and Brown correlation
% fluid = 4: multi-phase gas-oil-water flow, Mukherjee and Brill correlation
% fluid = 5: multi-phase gas-oil-water flow, Beggs and Brill correlation
% fluid = 6: multi-phase gas-oil-water flow, Shi et al. drift flux model (Matlab assignment)    
% 
% See script files 'example_flowline' and 'example_well' for examples of how to use this
% function.  
%
% alpha = inclination wrt. vertical, rad; alternatively alpha can be a survey file (matrix)
% with AHD values in the first column (in m) and inclination values in the second column (in
% rad).    
% d = inside diameter, m
% e = roughness, m
% oil = parameter to select black oil model or volatile oil table, -  
%   oil = 1: black oil; parameters computed with the aid of Standing correlations
%   oil = 2: black oil; parameters computed with the aid of Glaso correlations
%   oil = 3: volatile oil; parameters read from file 'vol_oil_table_01'
% p = [p_tot,p_grav,p_fric,p_acc], Pa
%   p_acc  = p_in + pressure increase due to acceleration losses, Pa
%   p_fric = p_in + pressure increase due to friction losses, Pa
%   p_grav = p_in + pressure increase due to head loss, Pa
%   p_tot  = p_in + pressure increase due to gravity, friction and acceleration losses, Pa
% p_in = pressure at s_in, Pa 
% p_out = pressure at s_out, Pa
% q_sc = [q_g_sc,q_o_sc,q_w_sc], m^3/s
%   q_g_sc = gas flow rate at standard conditions, m^3/s. Not relevant for fluid 1 
%   q_o_sc = oil flow rate at standard conditions, m^3/s. Not relevant for fluid 2 
%   q_w_sc = water flow rate at standard conditions, m^3/s. Not relevant for fluids 1 and 2 
%   Note: Flowrates in a production well need to have a negative value.  
% rho_sc = [rho_g_sc,rho_o_sc,rho_w_sc], kg/m^3 
%   rho_g_sc = gas density at standard conditions, kg/m^3. Not relevant for fluid 1
%   rho_o_sc = oil density at standard conditions, kg/m^3. Not relevant for fluid 2
%   rho_w_sc = water density at standard conditions, kg/m^3. Not relevant for fluids 1 and 2
% s = co-ordinate running from the separator to the reservoir, m
% s_in = starting point for the integration
% s_out = end point for the integration
% T_in = temperature at s_in, deg. C 
% T_out = temperature at s_out, deg. C 
%
% JDJ 11-03-02, last revised 11-04-16

interval = [s_in,s_out]; % integration interval, m
boundcon = [p_in,p_in,p_in,p_in]; % boundary condition, Pa 
options = []; % dummy variable, -
% options = odeset('MaxStep',10,'RelTol',1e-3); % tight tolerances; time consuming!
% options = odeset('MaxStep',1,'RelTol',1e-4); % more tight tolerances; more time consuming!
switch fluid
    case 1 % oil
        [s,p] = ode45('oil_dpds',interval,boundcon,options,alpha,d,e,oil,q_sc,rho_sc,...
            s_in,s_out,T_in,T_out);
    case 2 % gas
        [s,p] = ode45('gas_dpds',interval,boundcon,options,alpha,d,e,q_sc,rho_sc,...
            s_in,s_out,T_in,T_out);
    case 3 % multi-phase, Hagedorn and Brown
        [s,p] = ode45('Hag_Brown_dpds',interval,boundcon,options,alpha,d,e,oil,q_sc,...
            rho_sc,s_in,s_out,T_in,T_out);
    case 4 % multi-phase, Mukherjee and Brill
        [s,p] = ode45('Muk_Brill_dpds',interval,boundcon,options,alpha,d,e,oil,q_sc,...
            rho_sc,s_in,s_out,T_in,T_out);
    case 5 % multi-phase, Beggs and Brill
        [s,p] = ode45('Beggs_Brill_dpds',interval,boundcon,options,alpha,d,e,oil,q_sc,...
            rho_sc,s_in,s_out,T_in,T_out);
    case 6 % multi-phase, drift flux
        [s,p] = ode45('Drift_Flux_dpds',interval,boundcon,options,alpha,d,e,oil,q_sc,...
            rho_sc,s_in,s_out,T_in,T_out);
end
n = length(p);
p_out = p(n,1);
