function[x1] = impact_solver(t,x)
%function to solve newtons hypothesis alongside momentum conservation and
%constraint equations to return post impact state-space vector to feed into
%ode 89 again 
xi = x(length(t),1:end); % pre impact state-space in C2 basis
C2cO = Euler2DCM(xi(10:13)); %rotation matrix to get to C2 space from O space
OcC2 = transpose(C2cO); % rotation matrix to get to O space from C2 space
Icm_c2 = PInertiaTens(const.m,const.r0,const.h0);
I = OcC2*Icm_c2*C2cO;
psi = R(3,3); 
phi = impact_point(x(length(t),1:end)); %angle to define the polar co-ords of the point of impact
a11 =

A = []; %main part, use symbolic toolbox in a different script to generate the coeffs of the matrix
b = [];%main part, use sym toolbox in a diff script to get values
xo = linsolve(A,B); % moneyshot
%convert v and w into C2 space concat vector taking unchanged values and
%return function output



end

%%
