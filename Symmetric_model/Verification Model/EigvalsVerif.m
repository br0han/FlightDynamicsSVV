clear; clc;
run("Cit_par1.m");
%% Full Char Equation
A1 = 4*muc^2*KY2*(CZadot - 2*muc);
B1 = Cmadot*2*muc*(CZq + 2*muc) - Cmq*2*muc*(CZadot - 2*muc) - 2*muc*KY2*(CXu*(CZadot - 2*muc) - 2*muc*CZa);
C1 = Cma*2*muc*(CZq + 2*muc) - Cmadot*(2*muc*CX0 + CXu*(CZq + 2*muc)) + Cmq*(CXu*(CZadot - 2*muc) - 2*muc*CZa) + 2*muc*KY2*(CXa*CZu - CZa*CXu);
D1 = Cmu*(CXa*(CZq + 2*muc) - CZ0*(CZadot - 2*muc)) - Cma*(2*muc*CX0 + CXu*(CZq + 2*muc)) + Cmadot*(CX0*CXu - CZ0*CZu) + Cmq*(CXu*CZa - CZu*CXa);
E1 = -Cmu*(CX0*CXa + CZ0*CZa) + Cma*(CX0*CXu + CZ0*CZu);

lambda_full = roots([A1 B1 C1 D1 E1])
%% Simplified
A1 = 2*muc*KY2*(2*muc - CZadot);
B1 = -2*muc*KY2*CZa - (2*muc + CZq)*Cmadot - (2*muc - CZadot)*Cmq;
C1 = CZa*Cmq - (2*muc + CZq)*Cma;

lambda_simpSP = roots([A1, B1, C1])

A2 = 2*muc*(CZa*Cmq - 2*muc*Cma);
B2 = 2*muc*(CXu*Cma - Cmu*CXa) + Cmq*(CZu*CXa - CXu*CZa);
C2 = CZ0*(Cmu*CZa - CZu*Cma);

lambda_simpPhu = roots([A2, B2, C2])