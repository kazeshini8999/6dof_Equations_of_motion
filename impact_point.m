function [phi] = impact_point(x)
%returns the angle phi which denotes the position of impact on the
%coin,which is needed to compute the post impact state space
syms phi0
r_A_C2_C2 = [const.r0*cos(phi0);const.r0*sin(phi0);-const.h0/2];
r_C2_O_C2 = [x(1);x(2);x(3)];
r_A_O_C2  = r_A_C2_C2+r_C2_O_C2;
e = [x(10);x(11);x(12);x(13)];
r_A_O_O = transpose(Euler2DCM(e))*r_A_O_C2;
phi = double(solve(r_A_O_O(3),phi0)); % equates the z co-ord to zero in inertial space
end