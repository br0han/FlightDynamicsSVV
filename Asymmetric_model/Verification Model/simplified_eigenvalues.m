
run("ReaderParams.m");

% these parameters are for the dutch roll 
A = 8* mub^2 * KZ2;
B = -2* mub*(Cnr+ 2*KZ2 * CYb);
C = 4* mub * Cnb + CYb * Cnr;
% these parameters are for the dutch roll and the aperiodic roll 
A_da = 4*mub^2 *(KX2 * KZ2 - KXZ^2);
B_da = -mub *((Clr + Cnp)*KXZ + Cnr * KX2 + Clp *KZ2);
C_da = 2*mub*(Clb * KXZ + Cnb * KX2) + 1/4 * ( Clp * Cnr - Cnp * Clr);
D_da = 1/2 *(Clb * Cnp - Cnb * Clp);

lamda_dutch = roots([A B C])*(V0/b);
lamda_spiral = 2*CL*(Clb * Cnr - Cnb * Clr)/(Clp*( CYb * Cnr + 4*mub * Cnb) - Cnp*(CYb * Clr + 4*mub*Clb));
lamda_dutch_aperiodic = roots([A_da B_da C_da D_da]);

