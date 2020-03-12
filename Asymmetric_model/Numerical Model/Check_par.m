%hp0    = 5000;      	  % pressure altitude in the stationary flight condition [m]
V    = 59.9;            % true airspeed in the stationary flight condition [m/sec]
%m      = 9500;         % Kg

% rho0   = 1.2250;          % air density at sea level [kg/m^3] 
% lambda = -0.0065;         % temperature gradient in ISA [K/m]
% Temp0  = 288.15;          % temperature at sea level in ISA [K]
% R      = 287.05;          % specific gas constant [m^2/sec^2K]
% g      = 9.81;            % [m/sec^2] (gravity constant)
% 
% rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
% W      = m*g;
b = 13.36;
S = 24.2;
mu_b = 15.5;
C_L = 1.1360;
% V = V0;



Cy_bt_dot = 0;
Cy_bt = -0.9896;
Cy_p = -0.0870;
Cy_r = +0.430;
Cy_da = 0;
Cy_dr = +0.3037;

K_x = sqrt(0.012);
K_xz = 0.002;
K_z = sqrt(0.037);


Cn_bt_dot = 0;
Cn_bt = +0.1638;
Cn_p = -0.0108;
Cn_r = -0.1930;
Cn_da = -0.0286;
Cn_dr = -0.1261;

Cl_bt = -0.0772;
Cl_p = -0.3444;
Cl_r = 0.28;
Cl_da = -0.2349;
Cl_dr = 0.0286;