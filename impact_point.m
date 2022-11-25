function [phi] = impact_point(t,x)
%returns the angle phi which denotes the position of impact on the
%coin,which is needed to compute the post impact state space
syms phi0
xi = transpose(x(length(t),1:end));
e = [xi(10);xi(11);xi(12);xi(13)];
R = Euler2DCM(e);
cphi = R(3,3);
r_A_C2_C2 = [const.r0*cos(phi0);const.r0*sin(phi0);-sign(cphi)*const.h0/2];
r_C2_O_C2 = [xi(1);xi(2);xi(3)];
r_A_O_C2  = r_A_C2_C2+r_C2_O_C2;
r_A_O_O = transpose(Euler2DCM(e))*r_A_O_C2;
phi = double(solve(r_A_O_O(3),phi0)); % equates the z co-ord to zero in inertial space
end