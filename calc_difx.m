function [u v] = calc_difx(u0, v0, f0, fk, tstep)
global K_f0x K_Cfx K_u0x K_v0x Kv_f0x Kv_Cfx Kv_u0x Kv_v0x
Cf = (fk - f0) / tstep;
u = u0 * K_u0x + v0 * K_v0x + f0 * K_f0x + Cf * K_Cfx;
v = u0 * Kv_u0x + v0 * Kv_v0x + f0 * Kv_f0x + Cf * Kv_Cfx;