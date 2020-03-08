clear; close all; clc;
fd = importdata('C:\Users\Gebruiker\Documents\MATLAB\AE3212SA\matlab.mat');
t = fd.time.data;
%% Start time of eigen motions 
eigenMotions = ["Phugoid", "Aperiodic Roll", "Short Peroid", ...
				"Dutch Roll", "Dutch Roll Damped", "Sprial"];
aperi_roll_st = 59*60 + 10;					% Start time aperiodic roll [s]
short_peroid_st = 1*3600 + 0*60 + 35;		% Start time short period [s]
dutch_roll_st = 1*3600 + 1*60 + 57;			% Start time Dutch roll [s]
dutch_roll_damp_st = 1*3600 + 2*60 + 47;	% Start time Damped Dutch roll [s]
sprial_st = 1*3600 + 5*60 + 20;				% Start time sprial [s]
phugoid_st = 53*60 + 57;					% Start time phugoid [s]

motion_st = [phugoid_st, aperi_roll_st, short_peroid_st, ...
	dutch_roll_st, dutch_roll_damp_st, sprial_st];

t_idx = zeros(length(motion_st), 2);
buffer = [100 50 70 35 140 -170];
for i = 1:length(t_idx)
	idx = find(t >= motion_st(i));
	t_idx(i, 1) = idx(1);
	if i ~= length(t_idx)
		idx2 = find(t >= motion_st(i + 1) - buffer(i));
		t_idx(i, 2) = idx2(1);
	else
		idx2 = find(t >= motion_st(i) - buffer(i));
		t_idx(i, 2) = idx2(1);
	end
end
%% Flight Data
%airspeed h, alpha, theta, pitch rate, roll rate, yaw rate, loadfactor, 
TAS = fd.Dadc1_tas.data*0.51444444444;
CAS = fd.Dadc1_cas.data;
h = fd.Dadc1_alt.data;
AOA = fd.vane_AOA.data;
theta = fd.Ahrs1_Pitch.data;
phi = fd.Ahrs1_Roll.data;
pitch_rate = fd.Ahrs1_bPitchRate.data;
roll_rate = fd.Ahrs1_bRollRate.data;
yaw_rate = fd.Ahrs1_bYawRate.data;

%% Plots
% Phugoid
figure();
i = 1;
subplot(2, 1, 1)
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), TAS((t_idx(i, 1): t_idx(i, 2))));
ylabel("true airspeed [m/s]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), h((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("time [s]")
ylabel("altitude [m]")

subplot(2, 1, 2)
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), AOA((t_idx(i, 1): t_idx(i, 2))));
ylabel("angle [degrees]")
hold on
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), theta((t_idx(i, 1): t_idx(i, 2))));
grid on
axis tight
xlabel("time [s]")


% Aperi Roll
figure();
i = 2;
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), roll_rate((t_idx(i, 1): t_idx(i, 2))));
ylabel("Roll Rate [degree/s]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), yaw_rate((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("time [s]")
ylabel("Yaw Rate [degree/s]")

% Short Period
figure();
i = 3;
subplot(2, 1, 1)
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), TAS((t_idx(i, 1): t_idx(i, 2))));
ylabel("true airspeed [m/s]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), AOA((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("time [s]")
ylabel("Angle of Attack [degree/s]")

subplot(2, 1, 2)
yyaxis left
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), theta((t_idx(i, 1): t_idx(i, 2))));
ylabel("Theta [degrees]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), pitch_rate((t_idx(i, 1): t_idx(i, 2))));
grid on
axis tight
xlabel("time [s]")
ylabel("Pitch Rate [degrees/s]")

% Dutch Roll
figure();
i = 4;
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), roll_rate((t_idx(i, 1): t_idx(i, 2))));
ylabel("Roll Rate [degree/s]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), yaw_rate((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("time [s]")
ylabel("Yaw Rate [degree/s]")

% Damped Dutch Roll
figure();
i = 5;
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), roll_rate((t_idx(i, 1): t_idx(i, 2))));
ylabel("Roll Rate [degree/s]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), yaw_rate((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("time [s]")
ylabel("Yaw Rate [degree/s]")

% Sprial
figure();
i = 6;
subplot(2, 1, 1)
yyaxis left
title(eigenMotions(i))
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), TAS((t_idx(i, 1): t_idx(i, 2))));
ylabel("true airspeed [m/s]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), AOA((t_idx(i, 1): t_idx(i, 2))));
axis tight
grid on
xlabel("time [s]")
ylabel("Angle of Attack [degree/s]")

subplot(2, 1, 2)
yyaxis left
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), theta((t_idx(i, 1): t_idx(i, 2))));
ylabel("Theta [degrees]")

yyaxis right
plot(t(t_idx(i, 1): t_idx(i, 2)) - t(t_idx(i, 1)), phi((t_idx(i, 1): t_idx(i, 2))));
grid on
axis tight
xlabel("time [s]")
ylabel("Bank angle [degrees]")
