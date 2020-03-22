clear; clc; close all;
[FD, eigenMotions] = FlightData();			% Get Flight Data (FD) attributes
InitMass = (9165 + 4100)*0.45359237 + 760;	% Take off mass
%% Symmetric Eigenmotions
Symmetric = ["Phugoid", "Short Peroid"];	% Name of Symmetric Eigenmotions (case sensitive)
IDXs = ismember(eigenMotions, Symmetric);	% Slice columns associated to Symmetric motions
sym_eigenMotions = eigenMotions(IDXs);
t_idxs = FD.t_idx(IDXs, :);

for i = 1:length(sym_eigenMotions)
	% Compute Responce of State Space model from recorded stick input 
	idx0 = t_idxs(i, 1); idx1 = t_idxs(i, 2);
	input = (FD.delta_e(idx0: idx1) - FD.delta_e(idx0));
	mass = InitMass - FD.fuel_used_tot(idx0);
	init = [0, 0, 0, 0];
	[y, time] = SymmStateSpace(FD.h(idx0), FD.TAS(idx0), FD.AOA(idx0), FD.theta(idx0), mass, ...
							   FD.pitch_rate(idx0), FD.t(idx0: idx1) - FD.t(idx0), input, init);
	% Plots
	figure();
	subplot(4, 1, 1);
	plot(time, y(:, 1));
	title(Symmetric(i));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("u [m/s]", "Interpreter", "latex")
	hold on
	plot(time, FD.TAS(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on
	
	subplot(4, 1, 2);
	plot(time, y(:, 2));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\alpha$$ [rad]", "Interpreter", "latex")
	ylim([0, max(y(:, 2)) + 0.05])
	hold on
	plot(time, FD.AOA(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on

	subplot(4, 1, 3);
	plot(time, y(:, 3));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\theta$$ [rad]", "Interpreter", "latex")
	hold on
	plot(time, FD.theta(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on

	subplot(4, 1, 4);
	plot(time, y(:, 4));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\dot{\theta}$$ [rad/s]", "Interpreter", "latex");
	hold on
	plot(time, FD.pitch_rate(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on
	
	figure();
	plot(time, input)
	title([Symmetric(i), "Elevator Input"])
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\delta_{e}$$ [rad]", "Interpreter", "latex");
	grid on
end
%% Asymmetric Eigenmotions
Asymmetric = ["Aperiodic Roll", "Dutch Roll", "Spiral"];		% Name of Asymmetric Eigenmotions (case sensitive)
IDXa = ismember(eigenMotions, Asymmetric);						% Slice columns associated to Asymmetric motions
asym_eigenMotions = eigenMotions(IDXa);
t_idxa = FD.t_idx(IDXa, :);
delta_a0 = 0.01057;		delta_r0 = -0.00659;					% Zero Error in aileron and rudder measurement
for i = 1:length(asym_eigenMotions)
	% Compute Responce of State Space model from recorded aileron and rudder input
	idx0 = t_idxa(i, 1); idx1 = t_idxa(i, 2);
	input = [FD.delta_a(idx0: idx1) - delta_a0, FD.delta_r(idx0: idx1) - delta_r0];
	mass = InitMass - FD.fuel_used_tot(idx0);
	init = [0, 0, 0, 0];
	[y, time] = ASymmStateSpace(FD.h(idx0), FD.TAS(idx0), FD.AOA(idx0), FD.theta(idx0), 0, FD.phi(idx0), FD.roll_rate(idx0), ...
								FD.yaw_rate(idx0), mass, FD.t(idx0: idx1) - FD.t(idx0), input, init);
	% Plots
	figure();
	subplot(4, 1, 1);
	plot(time, y(:, 1));
	title(asym_eigenMotions(i));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\beta$$ [rad]", "Interpreter", "latex")
	legend("Simulation")
	grid on
	
	subplot(4, 1, 2);
	plot(time, y(:, 2));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\phi$$ [rad]", "Interpreter", "latex")
	hold on
	plot(time, FD.phi(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on

	subplot(4, 1, 3);
	plot(time, y(:, 3));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\dot{\phi}$$ [rad/s]", "Interpreter", "latex")
	hold on
	plot(time, FD.roll_rate(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on

	subplot(4, 1, 4);
	plot(time, y(:, 4));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\dot{\psi}$$ [rad/s]", "Interpreter", "latex")
	hold on
	plot(time, FD.yaw_rate(idx0: idx1));
	legend("Simulation", "Flight Data")
	grid on
	
	figure();
	subplot(2, 1, 1)
	plot(time, input(:, 1))
	title([Asymmetric(i), "Aileron Input"])
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\delta_{a}$$ [rad]", "Interpreter", "latex");
	grid on
	
	subplot(2, 1, 2)
	plot(time, input(:, 2))
	title([Asymmetric(i), "Rudder Input"])
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\delta_{r}$$ [rad]", "Interpreter", "latex");
	grid on
end