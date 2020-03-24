clear; clc; close all;
%% Start time of eigenmotions (Reference Data)
fd = importdata('matlab.mat');
aperi_roll_st = 59*60 + 10;							% Start time aperiodic roll [s]
short_peroid_st = 1*3600 + 0*60 + 33;					% Start time short period [s]
dutch_roll_st = 1*3600 + 1*60 + 57;					% Start time Dutch roll [s]
dutch_roll_damp_st = 1*3600 + 2*60 + 47;				% Start time Damped Dutch roll [s]
spiral_st = 1*3600 + 5*60 + 20;						% Start time sprial [s]
phugoid_st = 53*60 + 51;								% Start time phugoid [s]
buffer = [211 35 12 12 12 160];

%% Start time of eigenmotions (Our Filight data)
% fd = importdata('FTISxprt-20200306_flight3.mat');
% aperi_roll_st = 1*3600 + 1*60 + 15;						% Start time aperiodic roll [s]
% short_peroid_st = 56*60 + 50;							% Start time short period [s]
% dutch_roll_st = 57*60 + 94;								% Start time Dutch roll [s]
% dutch_roll_damp_st = 59*60 + 20;						% Start time Damped Dutch roll [s]
% spiral_st = 1*3600 + 3*60 + 53;							% Start time sprial [s]
% phugoid_st = 53*60 + 20;								% Start time phugoid [s]
% buffer = [175 35 20 15 15 180];

%% Index Slicing
t = fd.time.data;
eigenMotions = ["Phugoid", "Aperiodic Roll", "Short Peroid", ...
				"Dutch Roll", "Dutch Roll Damped", "Spiral"];
motion_st = [phugoid_st, aperi_roll_st, short_peroid_st, ...
	dutch_roll_st, dutch_roll_damp_st, spiral_st];

t_idx = zeros(length(motion_st), 2);
for i = 1:length(t_idx)
	idx = find(t >= motion_st(i));
	t_idx(i, 1) = idx(1);
	idx2 = find(t >= motion_st(i) + buffer(i));
	t_idx(i, 2) = idx2(1);
end
% Flight Data
% airspeed h, alpha, theta, pitch rate, roll rate, yaw rate
TAS = fd.Dadc1_tas.data*0.51444444444;
h = fd.Dadc1_alt.data*0.3048;
AOA = fd.vane_AOA.data*pi/180;
theta = fd.Ahrs1_Pitch.data*pi/180;
phi = fd.Ahrs1_Roll.data*pi/180;
pitch_rate = fd.Ahrs1_bPitchRate.data*pi/180;
roll_rate = fd.Ahrs1_bRollRate.data*pi/180;
yaw_rate = fd.Ahrs1_bYawRate.data*pi/180;
delta_e = fd.delta_e.data*pi/180;
delta_a = fd.delta_a.data*pi/180;

% DATA for eigenvalues
i = 1; %phugoid
Phugoid_V = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), TAS(t_idx(i, 1): t_idx(i, 2))];
Phugoid_th = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), theta(t_idx(i, 1): t_idx(i, 2))];
TAS_data_P = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), TAS(t_idx(i, 1): t_idx(i, 2))];

i = 4; %Dutch roll
Dutch_rollrate = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), roll_rate(t_idx(i, 1): t_idx(i, 2))];
Dutch_yawrate = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), yaw_rate(t_idx(i, 1): t_idx(i, 2))];
TAS_data_D = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), TAS(t_idx(i, 1): t_idx(i, 2))];


i = 5; %Damped Dutch roll
D_Dutch_rollrate = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), roll_rate(t_idx(i, 1): t_idx(i, 2))];
D_Dutch_yawrate = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), yaw_rate(t_idx(i, 1): t_idx(i, 2))];
TAS_data_DD = [transpose(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1))), TAS(t_idx(i, 1): t_idx(i, 2))];

%% Plots
% Phugoid
figure();
i = 1;
subplot(2, 1, 1)
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), TAS((t_idx(i, 1): t_idx(i, 2))));
ylabel("u [m/s]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), h((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("t [s]", "Interpreter", "latex")
ylabel("h [m]", "Interpreter", "latex")

subplot(2, 1, 2)
yyaxis left
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), AOA((t_idx(i, 1): t_idx(i, 2))));
ylabel("$$\alpha$$ [rad]", "Interpreter", "latex")
ylim([0, max(AOA((t_idx(i, 1): t_idx(i, 2)))) + 0.05])
yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), theta((t_idx(i, 1): t_idx(i, 2))));
grid on
ylabel("$$\theta$$ [rad]", "Interpreter", "latex")
axis tight
xlabel("t [s]", "Interpreter", "latex")


% Aperi Roll
figure();
i = 2;
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), roll_rate((t_idx(i, 1): t_idx(i, 2))));
ylabel("$$\dot{\phi}$$ [rad/s]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), yaw_rate((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("t [s]", "Interpreter", "latex")
ylabel("$$\dot{\psi}$$ [rad/s]", "Interpreter", "latex")

% Short Period
figure();
i = 3;
subplot(2, 1, 1)
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), TAS((t_idx(i, 1): t_idx(i, 2))));
ylabel("u [m/s]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), AOA((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("t [s]", "Interpreter", "latex")
ylabel("$$\alpha$$ [rad]", "Interpreter", "latex")

subplot(2, 1, 2)
yyaxis left
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), theta((t_idx(i, 1): t_idx(i, 2))));
ylabel("$$\theta$$ [rad]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), pitch_rate((t_idx(i, 1): t_idx(i, 2))));
grid on
axis tight
xlabel("t [s]", "Interpreter", "latex")
ylabel("$$\dot{\theta}$$ [rad/s]", "Interpreter", "latex")

% Dutch Roll
figure();
i = 4;
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), roll_rate((t_idx(i, 1): t_idx(i, 2))));
ylabel("$$\dot{\phi}$$ [rad/s]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), yaw_rate((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("t [s]", "Interpreter", "latex")
ylabel("$$\dot{\psi}$$ [rad/s]", "Interpreter", "latex")

% Damped Dutch Roll
figure();
i = 5;
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), roll_rate((t_idx(i, 1): t_idx(i, 2))));
ylabel("$$\dot{\phi}$$ [rad/s]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), yaw_rate((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("t [s]", "Interpreter", "latex")
ylabel("$$\dot{\psi}$$ [rad/s]", "Interpreter", "latex")

% Spiral
figure();
i = 6;
subplot(2, 1, 1)
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), TAS((t_idx(i, 1): t_idx(i, 2))));
ylabel("u [m/s]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), h((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("t [s]", "Interpreter", "latex")
ylabel("h [m]", "Interpreter", "latex")

subplot(2, 1, 2)
yyaxis left
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), theta((t_idx(i, 1): t_idx(i, 2))));
ylabel("$$\theta$$ [rad]", "Interpreter", "latex")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), phi((t_idx(i, 1): t_idx(i, 2))));
grid on
axis tight
xlabel("t [s]", "Interpreter", "latex")
ylabel("$$\psi$$ [rad]", "Interpreter", "latex")
