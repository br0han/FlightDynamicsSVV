run("ReaderParams.m");



A = 16*(mub^3)*(KX2*KZ2 - KXZ^2);

B = -4*(mub^2)*(2*CYb*(KX2*KZ2 - KXZ^2) + Cnr*KX2 + Clp*KZ2 + (Clr + Cnp)*(KXZ));

C = 2*mub*((CYb*Cnr - CYr*Cnb)*KX2 + (CYb*Clp - Clb*CYp)*KZ2 + ((CYb*Cnp - Cnb*CYp) + (CYb*Clr - Clb*CYr))*KXZ + 4*mub*Cnb*KX2 + 4*mub*Clb*KXZ + (1/2)*(Clp*Cnr - Cnp*Clr));

D = -4*mub*CL*(Clb*KZ2 + Cnb*KXZ) + 2*mub*(Clb*Cnp - Cnb*Clp) + (1/2)*(CYb*(Clr*Cnp - Cnr*Clp)) + (1/2)*(CYp*(Clb*Cnr - Cnb*Clr)) + (1/2)*(CYr*(Clp*Cnb - Cnp*Clb));

E = CL*(Clb*Cnr - Cnb*Clr);

R = B*C*D - A*(D^2) - (B^2)*E;
 
lambdabs = roots([A B C D E])*(V0/b)

