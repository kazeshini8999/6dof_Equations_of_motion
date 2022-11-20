function [I] = PInertiaTens(m,r,h)
%PInertiaTens returns the principle moment of inertia tensor for a cylinder in the C2
%frame basis, or the principle axes frame basis
m1 = (1/12*const.m*const.h0^2 + 1/4*const.m*const.r0^2);
m2 = m1;
m3 = (1/2*const.m*const.r0^2);
I = [m1 0 0;0 m2 0;0 0 m3];

end