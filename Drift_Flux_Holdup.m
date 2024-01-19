% Holdup calculation using Hagedown-Brown correlation
function  [H_g,H_l,lambda_l,lambda_g] = Drift_Flux_Holdup(alpha,d,q_g,q_l,rho_l,rho_g,sigma_gl)
% Input variables (that are considered as standard, from assignment 3)
% alpha = 60; % wellbore inclination , rad
% d = 0.0623;   % inside diameter, m
% e = 7.62e-5;  % roughness, m
% mu_g = 10e-6; % gas viscosity, Pa.s
% mu_l = 0.015; % liquid viscosity, Pa.s
% q_g_sc = -0.03; % Gas flowrate, m3/s
% q_o_sc = -0.01; % Oil flowrate, m3/s
% rho_o_sc = 762; % Oil density at standard conditions, kg/m3
% rho_g_sc = 69;  % Gas density at standard conditions, kg/m3
% sigma_gl = 0.01; % local gas-liquid interfacial tension, N/m 

% End of input data--------------------------------

% Computing superficial velocities
A = 3.1415*(d^2)/4; % Pipe area
v_sg = q_g/A;       % Superficial gas velocity, m/s
v_sl = q_l/A;       % Superficial liquid velocity, m/s
v_ms = v_sg + v_sl; % Superficial mixture velocity, m/s

% Computing dimensionless numbers:
g = 9.81;   % acceleration of gravity, m/s^2
N_d = d*sqrt(g*(rho_l-rho_g)/sigma_gl); % Dimensionless pipe diameter

    if N_d >= 2 && N_d <= 70
    K_u = (1.0152e-5*N_d^3 - 2.3396e-3*N_d^2 + 8.085e-1*N_d - 1.5934) / (1.9551e-1*N_d + 1); % Critical Kutateladze number
        else N_d < 2 || N_d > 70
            error('N_d is not in the range [2, 70]');
    end

% Critical and gas flooding velocities calculation
V_c = -(sigma_gl*g*(rho_l-rho_g)/(rho_l^2))^(1/4); % Critical velocity, m/s
V_g = K_u*sqrt(rho_l/rho_g)*V_c; % Gas flooding velocity, m/s

% Lambda g and lambda l
lambda_g = q_g/(q_l + q_g) % Gas no-slip holdup
lambda_l = 1 - lambda_g    % Liquid no-slip holdup

% Choosing optimal parameters for large and small pipe diameters:
    if d < 0.1  % for small pipe diameter
        c_o_bub = 1.2;
        beta_av = 0.6;
        a1 = 0.06;
        a2 = 0.12;
        m0 = 1.27;
        n1 = 0.24;
        n2 = 1.08;
    else % For large pipe diameter when d > 0.1 m
        c_o_bub = 1;
        beta_av = 1;
        a1 = 0.06;
        a2 = 0.21;
        m0 = 1.85;
        n1 = 0.21;
        n2 = 0.95;
    end

% Multiplier m_alpha to take account of wellbore inclination angle
m_alpha = m0*cos(alpha)^n1*(1+sin(alpha))^n2;

f_0 = 0.5; % Damping factor

% Calculating Holdup using Picard iterations
H_g_0 = lambda_g;
tolerance = 1e-4;

    for i = 1:1000
        beta = max(H_g_0,H_g_0*(v_ms/V_g));
        gamma = min(max(0,(beta-beta_av)/(1-beta_av)),1);
        C_0 = c_o_bub/(1+(c_o_bub-1)*gamma^2);
           
        if H_g_0 < a1 % Determine which K value should be used
            K = 1.53/C_0;
            elseif H_g_0 > a2
            K = K_u;
                else
                K = 1.53/C_0 + (H_g_0-a1)*(K_u - 1.53/C_0)/(a2 - a1);
           end
    
       v_d = m_alpha*(1 - H_g_0*C_0)*C_0*K*V_c/(H_g_0*C_0*sqrt(rho_g/rho_l)+1-H_g_0*C_0);
       H_g = (1-f_0)*H_g_0 + f_0*v_sg/(v_d + C_0*v_ms); % gas hold-up
    
        error = abs(H_g - H_g_0);
        if error > tolerance
            H_g_0 = H_g; % Reality check (not included in [1])
        else 
           H_g = H_g_0
            break
        end
    
    end
    H_l = 1-H_g; % gas hold-up, -
if q_g > 0 && H_l>lambda_l
       warning('downward flow, H_l>lambda_l');
     H_l=lambda_l;
     H_g=1-H_l;

end