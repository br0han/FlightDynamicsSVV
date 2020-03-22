function [responce, t_responce] = ASymmStateSpace(hp0, V0, alpha0, th0, beta0, phi0, p, r, m, time, input, init)
	run("Cit_par.m");
	% Rewritten form 
	% C1*x_dot + C2*x + C3*u = 0
	dim_cnst = b/V0;

	C_1 = [(CYbdot -(2*mub))*dim_cnst,		0,				0,						0;				...
					  0,			  -dim_cnst/2,			0,						0;				...
					  0,					0,		-4*mub*(KX2)*dim_cnst,   4*mub*KXZ*dim_cnst;	...
			    Cnbdot*dim_cnst,			0,		 4*mub*KXZ*dim_cnst,	-4*mub*(KZ2)*dim_cnst];

	C_2 = [CYb,		CL,		CYp,    (CYr - 4*mub);	...
			0,		0,		 1,			  0;		...
		   Clb,     0,      Clp,		 Clr;		...
		   Cnb,		0,		Cnp,		 Cnr];

	C_3 = [CYda,	CYdr; ...
			0,		 0;   ...
		   Clda,	Cldr; ...
		   Cnda,	Cndr];

	% State Space System
	% y = x (chosen output vector is same as state vector)
	% Aa = -(C1^-1)*C2
	% Ba = - (C1^-1)*C1
	% Ca = Identity(4)
	% Da = null vector

	Aa = inv(C_1)*(-C_2);
	Ba = inv(C_1)*(-C_3);
	Ca = eye(4);
	Ca(3, 3) = 2/dim_cnst;
	Ca(4, 4) = 2/dim_cnst;
	Ca = Ca*-1;
	Da = zeros(4, 2);

	sys_a = ss(Aa, Ba, Ca, Da);
	sys_a.StateName = {'\beta', '\phi', 'p', 'r'};
	sys_a.InputName = {'Aileron', 'Rudder'};
	sys_a.OutputName = {'\beta', '\phi', 'p', 'r'};
	[responce, t_responce] = lsim(sys_a, input, time, init);
	responce = responce + [beta0, phi0, p, r];
end