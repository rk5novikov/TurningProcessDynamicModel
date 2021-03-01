function [Fz,Fr]=calc_force(Kc,s,h_z,h_r)
% Cutting forces calculation
Fz=0;
Fr=0;
for j=1:length(s)
    Fz=Fz+Kc*s(j)*h_z(j);
    Fr=Fr+Kc*s(j)*h_r(j); 
end