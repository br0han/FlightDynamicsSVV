%Asymmetric equations of Motion in State Space System

% Parameters

mu_b = 100;
C_L = 1;
b = 15.911;
V = 100;
S = 30;


Cy_bt_dot = 0;
Cy_bt = -0.7500;
Cy_p = -0.0304;
Cy_r = +0.8495;
Cy_da = -0.0400;
Cy_dr = +0.2300;

K_x = sqrt(0.019);
K_xz = 0.002;
K_z = sqrt(0.042);


Cn_bt_dot = 0;
Cn_bt = +0.1348;
Cn_p = -0.0602;
Cn_r = -0.2061;
Cn_da = -0.0120;
Cn_dr = -0.0939;

Cl_bt = -0.10260;
Cl_p = -0.71085;
Cl_r = +0.23760;
Cl_da = -0.23088;
Cl_dr = +0.03440;


% Rewritten form 
% C1*x_dot + C2*x + C3*u = 0

C_1 = [2*(Cy_bt_dot - 2*mu_b),0,0,0; 0,-1,0,0; 0,0,(-4*mu_b*(K_x^2)*(b/V)),(4*mu_b*K_xz*(b/V));
    2*Cn_bt_dot,0,(4*mu_b*K_xz*(b/V)),(-4*mu_b*(K_z^2)*(b/V))];

C_2 = [(Cy_bt*(2*V/b)),(C_L*(2*V/b)),Cy_p,((Cy_r - (4*mu_b))); 0,0,1,0; Cl_bt*(2*V/b),0,Cl_p,Cl_r;
    (Cn_bt*(2*V/b)),0,Cn_p,Cn_r];

C_3 = [-Cy_da*(2*V/b), -Cy_dr*(2*V/b); 0,0 ; -Cl_da*(2*V/b), -Cl_dr*(2*V/b) ; -Cn_da*(2*V/b), -Cn_dr*(2*V/b)];


% State Space System

% y = x (chosen output vector is same as state vector)
% Aa = -(C1^-1)*C2
% Ba = - (C1^-1)*C3
% Ca = Identity(4)
% Da = null vector

Aa = -inv(C_1)*C_2;
Ba = -inv(C_1)*C_3;
Ca = eye(4);
Da = 0;

sys_a = ss(Aa,Ba,Ca,Da);
step(sys_a);

