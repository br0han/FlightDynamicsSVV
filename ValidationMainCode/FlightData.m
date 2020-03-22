function [FD, eigenMotions] = FlightData()
	fd = importdata('FTISxprt-20200306_flight3.mat');
	t = fd.time.data;
	% Start time of eigen motions 
	eigenMotions = ["Phugoid", "Aperiodic Roll", "Short Peroid", ...
					"Dutch Roll", "Dutch Roll Damped", "Spiral"];

	aperi_roll_st = 1*3600 + 1*60 + 15;						% Start time aperiodic roll [s]
	short_peroid_st = 56*60 + 55;							% Start time short period [s] p=44.0333333333 s
	dutch_roll_st = 57*60 + 94;								% Start time Dutch roll [s]
	dutch_roll_damp_st = 59*60 + 20;						% Start time Damped Dutch roll [s]
	spiral_st = 1*3600 + 3*60 + 53;							% Start time sprial [s]
	phugoid_st = 53*60 + 25;								% Start time phugoid [s]
% 	phugoid_st = 52*60;
	motion_st = [phugoid_st, aperi_roll_st, short_peroid_st, ...
		dutch_roll_st, dutch_roll_damp_st, spiral_st];

	t_idx = zeros(length(motion_st), 2);

	buffer = [175 25 20 15 15 80];

	for i = 1:length(t_idx)
		idx = find(t >= motion_st(i));
		t_idx(i, 1) = idx(1);
		idx2 = find(t >= motion_st(i) + buffer(i));
		t_idx(i, 2) = idx2(1);
	end
	
	% Flight Data Attributes
	FD.t = t;
	FD.t_idx = t_idx;
	FD.TAS = fd.Dadc1_tas.data*0.51444444444;
	FD.h = fd.Dadc1_alt.data*0.3048;
	FD.AOA = fd.vane_AOA.data*pi/180;
	FD.theta = fd.Ahrs1_Pitch.data*pi/180;
	FD.phi = fd.Ahrs1_Roll.data*pi/180;
	FD.pitch_rate = fd.Ahrs1_bPitchRate.data*pi/180;
	FD.roll_rate = fd.Ahrs1_bRollRate.data*pi/180;
	FD.yaw_rate = fd.Ahrs1_bYawRate.data*pi/180;
	FD.delta_e = fd.delta_e.data*pi/180;
	FD.delta_a = fd.delta_a.data*pi/180;
	FD.delta_r = fd.delta_r.data*pi/180;
	FD.fuel_used_tot = (fd.lh_engine_FU.data + fd.rh_engine_FU.data)*0.453592;
end