hp0    = 6900;      	  % pressure altitude in the stationary flight condition [m]
V0    = 145;            % true airspeed in the stationary flight condition [m/sec]
m      = 59020;         % Kg

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;

b = 37.49;
S = 153.5;
mub = 17.219;
CL = 2*W/(rho*(V0^2)*S);
% V = V0;



CYbdot = 0;
CYb = -0.596;
CYp = 0;
CYr = 0.369;
CYda = 0;
CYdr = 0.215;

KX2 = 0.0283;
KXZ = 0;
KZ2 = 0.0471;


Cnbdot = 0;
Cnb = 0.1173;
Cnp = -0.021;
Cnr = -0.18;
Cnda = 0.0052;
Cndr = -0.1030;

Clb = -0.1374;
Clp = -0.52;
Clr = 0.144;
Clda = -0.0975;
Cldr = 0;