
%run("ReaderParams.m");
run("Cit_par.m");
% these parameters are for the dutch roll 
A = 8* mub^2 * KZ2;
B = -2* mub*(Cnr+ 2*KZ2 * CYb);
C = 4* mub * Cnb + CYb * Cnr;
% these parameters are for the dutch roll and the aperiodic roll 
A_da = 4*mub^2 *(KX2 * KZ2 - KXZ^2);
B_da = -mub *((Clr + Cnp)*KXZ + Cnr * KX2 + Clp *KZ2);
C_da = 2*mub*(Clb * KXZ + Cnb * KX2) + 1/4 * ( Clp * Cnr - Cnp * Clr);
D_da = 1/2 *(Clb * Cnp - Cnb * Clp);

A_dutch2= -2*mub*KZ2;
B_dutch2= 0.5*Cnr;
C_dutch2 = -Cnb;

lamda_aperiodic_roll = Clp/(4*mub*KX2);
lamda_dutch = roots([A B C]);
lamda_spiral = 2*CL*(Clb * Cnr - Cnb * Clr)/(Clp*( CYb * Cnr + 4*mub * Cnb) - Cnp*(CYb * Clr + 4*mub*Clb));
lamda_dutch_aperiodic = roots([A_da B_da C_da D_da]);

lamda_dutch2 = roots([A_dutch2 B_dutch2 C_dutch2]);

% eigenvalues with no simplifications 



A1 = 4*muc^2*KY2*(CZadot - 2*muc);
B1 = Cmadot*2*muc*(CZq + 2*muc) - Cmq*2*muc*(CZadot - 2*muc) - 2*muc*KY2*(CXu*(CZadot - 2*muc) - 2*muc*CZa);
C1 = Cma*2*muc*(CZq + 2*muc) - Cmadot*(2*muc*CX0 + CXu*(CZq + 2*muc)) + Cmq*(CXu*(CZadot - 2*muc) - 2*muc*CZa) + 2*muc*KY2*(CXa*CZu - CZa*CXu);
D1 = Cmu*(CXa*(CZq + 2*muc) - CZ0*(CZadot - 2*muc)) - Cma*(2*muc*CX0 + CXu*(CZq + 2*muc)) + Cmadot*(CX0*CXu - CZ0*CZu) + Cmq*(CXu*CZa - CZu*CXa);
E1 = -Cmu*(CX0*CXa + CZ0*CZa) + Cma*(CX0*CXu + CZ0*CZu);

lambda = roots([A1 B1 C1 D1 E1]); 
