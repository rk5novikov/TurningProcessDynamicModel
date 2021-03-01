function calc_const(c1, p0_21, n1, p1, c2, p0_22, n2, p2, tstep)
global K_f0x K_Cfx K_u0x K_v0x Kv_f0x Kv_Cfx Kv_u0x Kv_v0x
global K_f0y K_Cfy K_u0y K_v0y Kv_f0y Kv_Cfy Kv_u0y Kv_v0y
t = tstep;
K_f0x = (-n1/p1*exp(-n1*t)*sin(p1*t)-exp(-n1*t)*cos(p1*t) + 1.0)/c1;
K_Cfx = (-2*n1/p0_21-(p1^2-n1^2)/p0_21/p1*exp(-n1*t)*sin(p1*t)+2*n1/p0_21*exp(-n1*t)*cos(p1*t) + t)/c1;
K_u0x = (exp(-n1*t)*(cos(p1*t)+n1/p1*sin(p1*t)));
K_v0x = (1.0/p1*exp(-n1*t)*sin(p1*t));
%-------------------------------------------------------
Kv_f0x = (exp(-n1*t)*sin(p1*t)*(p0_21/p1))/c1;
Kv_Cfx = (exp(-n1*t)*sin(p1*t)*(-n1/p1)+exp(-n1*t)*cos(p1*t)*(-1)+1)/c1;
Kv_u0x = -n1*K_u0x + exp(-n1*t)*(-p1*sin(p1*t)+n1*cos(p1*t));
Kv_v0x = -n1*K_v0x + exp(-n1*t)*cos(p1*t);

K_f0y = (-n2/p2*exp(-n2*t)*sin(p2*t)-exp(-n2*t)*cos(p2*t) + 1.0)/c2;
K_Cfy = (-2*n2/p0_22-(p2^2-n2^2)/p0_22/p2*exp(-n2*t)*sin(p2*t)+2*n2/p0_22*exp(-n2*t)*cos(p2*t) + t)/c2;
K_u0y = (exp(-n2*t)*(cos(p2*t)+n2/p2*sin(p2*t)));
K_v0y = (1.0/p2*exp(-n2*t)*sin(p2*t));
%-------------------------------------------------------
Kv_f0y = (exp(-n2*t)*sin(p2*t)*(p0_22/p2))/c2;
Kv_Cfy = (exp(-n2*t)*sin(p2*t)*(-n2/p2)+exp(-n2*t)*cos(p2*t)*(-1)+1)/c2;
Kv_u0y = -n2*K_u0y + exp(-n2*t)*(-p2*sin(p2*t)+n2*cos(p2*t));
Kv_v0y = -n2*K_v0y + exp(-n2*t)*cos(p2*t);