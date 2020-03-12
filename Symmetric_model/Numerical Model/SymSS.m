clear; close all; clc;
%% Citation Data
% run("Cit_par.m");				
%% Data from FD reader page 133 Table 5-1
alpha0 = 0; th0 = 0;
V0 = 59.9; m = 4547.8; c = 2.022;
S = 24.2;
lh = 5.5;
muc = 102.7;
KY2 = 0.980;
CX0 = 0;		CZ0 = -1.1360;		
CXu = -0.2199;	CZu = -2.2720;		Cmu = 0;
CXa = 0.4653;	CZa = -5.1600;		Cma = -0.4300;
CXadot = 0;		CZadot = -1.4300;	Cmadot = -3.7000;
CXq = 0;		CZq = -3.8600;		Cmq = -7.0400;
CXde = 0;		CZde = -0.6238;		Cmde = -1.5530;
%% State Space
Q = c/V0;
C_1 = [-2*muc*Q,		0,				0,		0; ...
		  0,	(CZadot - 2*muc)*Q,		0,		0; ...
		  0,			0,				-Q,		0; ...
		  0,		Cmadot*Q,			0,  -2*muc*KY2*Q];
	  
C_2 = [-CXu,	   -CXa,		-CZ0,			0;		...
	   -CZu,	   -CZa,		CX0,	-(CZq + 2*muc); ...
	     0,		    0,		     0,			   -1;		...
	   -Cmu,	   -Cma,		 0,		      -Cmq];
	
C_3 = [-CXde;	...
	   -CZde;	...
		0;		...
	   -Cmde];

As = inv(C_1)*C_2;
Bs = inv(C_1)*C_3;

Cs = [V0, 0, 0, 0; ...
	   0, 1, 0, 0; ...
	   0, 0, 1, 0; ...
	   0, 0, 0, 1/Q];
Ds = zeros(4, 1);
sys_s = ss(As, Bs, Cs, Ds, 'StateName', {'u' 'AOA' 'theta' 'pitch rate'}, 'InputName', 'Delta elevator');
%% Solution
init = [0, alpha0, th0, 0];			% Initial Condition
t = 0:0.1:150;						% Time Domain
step = -0.005*ones(1, length(t));	% Step input -0.005 rad of delta_e
[y, t] = lsim(sys_s, step, t, init);

subplot(4, 1, 1);
plot(t, y(:, 1) + V0);
xlabel("time [s]")
ylabel("u [m/s]")
grid on
subplot(4, 1, 2);
plot(t, y(:, 2));
xlabel("time [s]")
ylabel("AoA [rad]")
grid on

subplot(4, 1, 3);
plot(t, y(:, 3));
xlabel("time [s]")
ylabel("theta [rad]")
grid on

subplot(4, 1, 4);
plot(t, y(:, 4));
xlabel("time [s]")
ylabel("q [rad/s]")
grid on

lambda = eig(As)*Q