clear all
clc
%Asymmetric equations of Motion in State Space System

% Parameters
run("ReaderParams.m")

% Rewritten form 
% C1*x_dot + C2*x + C3*u = 0
dim_cnst = b/V0;

C_1 = [(CYbdot -(2*mub))*dim_cnst   ,  0    ,  0    ,  0    ;   ...
        0   ,  -dim_cnst/2  ,    0  ,  0    ; ...
        0   ,  0    ,  -4*mub*(KX2)*dim_cnst   ,   4*mub*KXZ*dim_cnst  ; ...
        Cnbdot*dim_cnst  , 0 ,  4*mub*KXZ*dim_cnst   , -4*mub*(KZ2)*dim_cnst];
    
C_2 = [CYb    ,   CL     ,    CYp    ,    (CYr - 4*mub) ;...
        0   ,  0    ,    1     ,   0   ;...
        Clb   ,    0    ,     Clp     ,   Clr     ;...
        Cnb   ,   0   ,   Cnp ,   Cnr     ];

C_3 = [ CYda ,     CYdr   ; ...
        0   ,   0   ;   ...
        Clda   ,   Cldr   ;   ...
        Cnda   ,   Cndr       ];



% State Space System

% y = x (chosen output vector is same as state vector)
% Aa = -(C1^-1)*C2
% Ba = - (C1^-1)*C1
% Ca = Identity(4)
% Da = null vector

Aa = inv(C_1)*(-C_2);
Ba = inv(C_1)*(-C_3);
Ca = eye(4);
Ca(3,3) = 2/dim_cnst;
Ca(4,4) = 2/dim_cnst;
Da = zeros(4,2);

step_size = 0.01;
tot_time = 20;

aileron_pulsetime = 1;
rudder_pulsetime = 1;

t = 0:step_size:tot_time;
x_0 = [0,0,0,0];
u = zeros(2,length(t));
% u(2,1:rudder_pulsetime/step_size) = 0.025;
u(1,1:10) = 1000;


sys_a = ss(Aa,Ba,Ca,Da);
sys_a.StateName = {'\beta','\phi','p','r'};
sys_a.InputName = {'Aileron','Rudder'};
sys_a.OutputName = {'\beta','\phi','p','r'};

% Eigenvalues
lambda = eig(Aa)

figure();
plot(lambda,'x')
xlabel('\xi')
ylabel('j \eta')

%Dynamic Response 
figure();
% lsim(sys_a,u,t);
impulse(sys_a,t)



