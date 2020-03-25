clear; close all; clc;
%% Initial Parameters
c = 2.0569;
b = 15.911;
hp0 = 1000;
V0 = 100;
alpha0 = 0.1;
beta0 = 0;
th0 = 0.15;
phi0 = 0;
m = 5000;
p0 = 0;
q0 = 0;
r0 = 0;
t = 0:0.1:60;
t_a = 0:0.1:20;
folder = "C:\Users\Gebruiker\Desktop\Uni\SVV\FlightDynamicsPlots\InitialValueProblem";
%% Symmetric Deviation
name_s = ["u", "alpha", "theta", "q"];
Ysym = struct(name_s(1), SymmStateSpace(hp0, V0, alpha0, th0, m, q0, t, zeros(1, length(t)), [10/V0, 0, 0, 0], "NonAdjusted"), ...
		  name_s(2), SymmStateSpace(hp0, V0, alpha0, th0, m, q0, t, zeros(1, length(t)), [0, 0.1, 0, 0], "NonAdjusted"),  ...
		  name_s(3), SymmStateSpace(hp0, V0, alpha0, th0, m, q0, t, zeros(1, length(t)), [0, 0, 0.1, 0], "NonAdjusted"),  ...
		  name_s(4), SymmStateSpace(hp0, V0, alpha0, th0, m, q0, t, zeros(1, length(t)), [0, 0, 0, 0.1*c/V0], "NonAdjusted")); 

for i = 1:length(name_s)
	responce = Ysym.(name_s(i));
	fig = figure();
	ax(1) = subplot(2, 2, 1);
	plot(t, responce(:, 1));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("TAS [m/s]", "Interpreter", "latex")
	grid on

	ax(2) = subplot(2, 2, 2);
	plot(t, responce(:, 2));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\alpha$$ [rad]", "Interpreter", "latex")
	grid on

	ax(3) = subplot(2, 2, 3);
	plot(t, responce(:, 3));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\theta$$ [rad]", "Interpreter", "latex")
	grid on

	ax(4) = subplot(2, 2, 4);
	plot(t, responce(:, 4));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\dot{\theta}$$ [rad/s]", "Interpreter", "latex");
	grid on
	ax(i).XColor = "red";
	ax(i).YColor = "red";
% 	saveas(fig, fullfile(folder, join(["InitialDisturbance", name_s(i)], "_")), "png")
end
%% Asymmetric Deviations
name_a = ["beta", "phi", "p", "r"];
Yasym = struct(name_a(1), ASymmStateSpace(hp0, V0, alpha0, th0, beta0, phi0, p0, r0, m, t_a, zeros(2, length(t_a)), [-0.1, 0, 0, 0], "NonAdjusted"), ...
			  name_a(2), ASymmStateSpace(hp0, V0, alpha0, th0, beta0, phi0, p0, r0, m, t_a, zeros(2, length(t_a)), [0, -0.1, 0, 0], "NonAdjusted"),  ...
			  name_a(3), ASymmStateSpace(hp0, V0, alpha0, th0, beta0, phi0, p0, r0, m, t_a, zeros(2, length(t_a)), [0, 0, -0.1*b/(2*V0), 0], "NonAdjusted"),  ...
			  name_a(4), ASymmStateSpace(hp0, V0, alpha0, th0, beta0, phi0, p0, r0, m, t_a, zeros(2, length(t_a)), [0, 0, 0, -0.1*b/(2*V0)], "NonAdjusted")); 
for i = 1:length(name_a)
	responce = Yasym.(name_a(i));
	fig = figure();
	ax(1) = subplot(2, 2, 1);
	plot(t_a, responce(:, 1));
	
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\beta$$ [rad]", "Interpreter", "latex")
	grid on

	ax(2) = subplot(2, 2, 2);
	plot(t_a, responce(:, 2));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\phi$$ [rad]", "Interpreter", "latex")
	grid on

	ax(3) = subplot(2, 2, 3);
	plot(t_a, responce(:, 3));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\dot{\phi}$$ [rad/s]", "Interpreter", "latex")
	grid on

	ax(4) = subplot(2, 2, 4);
	plot(t_a, responce(:, 4));
	xlabel("t [s]", "Interpreter", "latex")
	ylabel("$$\dot{\psi}$$ [rad/s]", "Interpreter", "latex");
	grid on
	ax(i).XColor = "red";
	ax(i).YColor = "red";
% 	saveas(fig, fullfile(folder, join(["InitialDisturbance", name_a(i)], "_")), "png")
end

