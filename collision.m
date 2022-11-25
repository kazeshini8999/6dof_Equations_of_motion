function [check,isterminal,direction] = collision(t,x)
%event function to detect impact with the floor and stop the ode solver, 
%is called each time stpe and executes when the check value is zero 
%
e = x(10:13); % Euler parameters at that time step
R1 = Euler2DCM(e); % rotation matrix to get to C2 frame from O frame
OcC2 = transpose(R1); %rotation matrix to get to O frame from C2 
phi = acos(R1(3,3)); % angle between the k axes of the two frames
posC2 = x(1:3); % position vector of com in C2 frame
posO = OcC2*posC2;  % position vector of com in O frame

check = posO(3) - (const.r0*sin(phi)+const.h0/2*cos(phi)); % position of coin center accounting for radius of coin
isterminal = 1; % Do I want to stop the code when event happens, 1 for yes 0 for no
direction = -1; % do I care from which direction the check hits 0, value of 0 means I dont care
end