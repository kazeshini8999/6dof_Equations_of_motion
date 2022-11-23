function[x1] = impact_solver(t,x)
%function to solve newtons hypothesis alongside momentum conservation and
%constraint equations to return post impact state-space vector to feed into
%ode 89 again 
xi = x(length(t),1:end);
R = Euler2DCM(xi(10:13));
OcC2 = transpose(R); % rotation matrix to get to O space from C2 space
psi = R(3,3);
phi = impact_point(x(length(t),1:end)); %angle to define the polar co-ords of the point of impact
va0C2 = [xi(4:6)] + vec2cross([xi(7:9)])*[const.r0*cos(phi);const.r0*sin(phi);-1*Sign(psi)*const.h0/2];
va0O = OcC2*va0C2; % vel of point A wrt to O in O basis (known)




end

%%
%solve for unknowns in inertial space, convert all the vectors to C2 space
%and output them as a function output