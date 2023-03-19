[ param ]

Vt = 0.1
CLd = 0.003
CL = 0.001
kon = 0.091
koff = 0.001
kin = 0.11
kout = 0.009
keRL = 0.003
R0 = 12
Vc = 0.05;

[ cmt ] L R RL T

[ pk ]

R_0 = R0;

double keL = CL/Vc;

F_L = 1/Vc;

[ des ]

dxdt_L = - kon * L * R + koff * RL - keL * L - (CLd/Vc) * L + (CLd/Vt) * T;

dxdt_R = kin - kout * R - kon * L * R + koff * RL;

dxdt_RL = kon * L * R - (koff + keRL) * RL;

dxdt_T = (CLd/Vc) * L - (CLd/Vt) * T;

$CAPTURE CP = L
