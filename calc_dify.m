function [u v] = calc_dify(u0, v0, f0, fk, tstep)
global K_f0y K_Cfy K_u0y K_v0y Kv_f0y Kv_Cfy Kv_u0y Kv_v0y
Cf = (fk - f0) / tstep;
u = u0 * K_u0y + v0 * K_v0y + f0 * K_f0y + Cf * K_Cfy;
v = u0 * Kv_u0y + v0 * Kv_v0y + f0 * Kv_f0y + Cf * Kv_Cfy;