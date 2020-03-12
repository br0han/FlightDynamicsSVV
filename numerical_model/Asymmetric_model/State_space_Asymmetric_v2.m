clear all
clc
%Asymmetric equations of Motion in State Space System
run('Check_par.m')

% Parameters

% hp0    = 5000;      	  % pressure altitude in the stationary flight condition [m]
% V0     = 100;            % true airspeed in the stationary flight condition [m/sec]
% m      = 9500;         % Kg
% 
% rho0   = 1.2250;          % air density at sea level [kg/m^3] 
% lambda = -0.0065;         % temperature gradient in ISA [K/m]
% Temp0  = 288.15;          % temperature at sea level in ISA [K]
% R      = 287.05;          % specific gas constant [m^2/sec^2K]
% g      = 9.81;            % [m/sec^2] (gravity constant)
% 
% rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
% W      = m*g;
% b = 15.911;
% S = 30;
% mu_b = m/(rho*S*b);
% C_L = 2*W/(rho*V0^2*S);
% V = V0;
% 
% 
% 
% Cy_bt_dot = 0;
% Cy_bt = -0.7500;
% Cy_p = -0.0304;
% Cy_r = +0.8495;
% Cy_da = -0.0400;
% Cy_dr = +0.2300;
% 
% K_x = sqrt(0.019);
% K_xz = 0.002;
% K_z = sqrt(0.042);
% 
% 
% Cn_bt_dot = 0;
% Cn_bt = +0.1348;
% Cn_p = -0.0602;
% Cn_r = -0.2061;
% Cn_da = -0.0120;
% Cn_dr = -0.0939;
% 
% Cl_bt = -0.10260;
% Cl_p = -0.71085;
% Cl_r = +0.23760;
% Cl_da = -0.23088;
% Cl_dr = +0.03440;
% 
% %Defining Constants

y_bt = (V/b)*(Cy_bt/(2*mu_b));
y_phi = (V/b)*(C_L/(2*mu_b));
y_p = (V/b)*(Cy_p/(2*mu_b));
y_r = (V/b)*((Cy_r -(4*mu_b))/(2*mu_b));
y_da = (V/b)*(Cy_da/(2*mu_b));
y_dr = (V/b)*(Cy_dr/(2*mu_b));

const_intm = 4*mu_b*(((K_x^2)*(K_z^2)) - (K_xz^2));
l_bt = (V/b)*((Cl_bt*(K_z^2)) + ((Cn_bt)*K_xz))/const_intm;
l_p = (V/b)*((Cl_p*(K_z^2)) + ((Cn_p)*K_xz))/const_intm;
l_r = (V/b)*((Cl_r*(K_z^2)) + ((Cn_r)*K_xz))/const_intm;
l_da = (V/b)*((Cl_da*(K_z^2)) + ((Cn_da)*K_xz))/const_intm;
l_dr = (V/b)*((Cl_dr*(K_z^2)) + ((Cn_dr)*K_xz))/const_intm;

n_bt = (V/b)*((Cl_bt*(K_xz)) + ((Cn_bt)*(K_x^2)))/const_intm;
n_p = (V/b)*((Cl_p*(K_xz)) + ((Cn_p)*(K_x^2)))/const_intm;
n_r = (V/b)*((Cl_r*(K_xz)) + ((Cn_r)*(K_x^2)))/const_intm;
n_da = (V/b)*((Cl_da*(K_xz)) + ((Cn_da)*(K_x^2)))/const_intm;
n_dr = (V/b)*((Cl_dr*(K_xz)) + ((Cn_dr)*(K_x^2)))/const_intm;

% State Space Matrix
A_as = [y_bt,y_phi,y_p,y_r; 0,0,(2*V/b),0; l_bt,0,l_p,l_r; n_bt,0,n_p,n_r];

B_as = [0,y_dr; 0,0; l_da,l_dr; n_da,n_dr];

C_as = eye(4);
C_as(3,3) = (2*V/b);
C_as(4,4) = (2*V/b);

D_as = zeros(4,2);

x_0 = [0,0,0,0];
t = linspace(0,50,100);
u = zeros(2,length(t));
u(1,:) = 1;

sys_as = ss(A_as,B_as,C_as,D_as);
sys_as.StateName = {'\beta','\phi','p','r'};
sys_as.InputName = {'Aileron','Rudder'};
sys_as.OutputName = {'\beta','\phi','p','r'};


figure()
lsim(sys_as,u,t,x_0);

eig(A_as)
