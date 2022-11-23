function [v2] = OcC2(v1,x) % you need to include t when you input x
%function to convert vectors from C2 frame basis to O frame basis
OcC2 = transpose(Euler2DCM(x(10:13)));
v2 = OcC2*v1;
end
