function [v2] = C2cO(v1,x) % you need to include t when you input x
%function to convert vectors from O frame basis to C frame basis
C2cO = Euler2DCM(x(10:13));
v2 = OcC2*v1;
end
