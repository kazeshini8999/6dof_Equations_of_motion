%% code for ode function calls
clc

% profile on
%% intial conditions
tspan = [0 0.1]; % event collision will detect when to stop integration, higher values lead to smoother curves
E0 = DCM2Euler(7*pi/180,0,0); %intial euler parameter matrix
r1 = (Euler2DCM(E0))*[0;0;0.41];
x0 = post2;
%x0 = [r1(1);r1(2);r1(3);0;0;0;40.15;0;0;E0(1);E0(2);E0(3);E0(4)]; %intial conditions
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-9,'Events',@collision); %ode settings
%% ode 45 call
[t,x]=ode45(@(t,x) dof(t,x),tspan,x0,options); %function call
%% system of ODES defination 
function xdot = dof(t,x)
% function to get system of odes
v0 = x(4:6); %velocity of com
w0 = x(7:9); % angular velocity of frame C2 wrt O
I = PInertiaTens(const.m,const.r0,const.h0); % principle inertia tensor matrix
Ixx = I(1,1); %principle moments
Iyy = I(2,2);
Izz = I(3,3);
e = x(10:13); %Euler parameters 
R1 = Euler2DCM(e); %Rotation matrix to get to C2 frame basis from O

f = drag1(v0,w0,const.r0,const.h0,const.Lt,const.Ln,const.b1)+R1*[0;0;-9.81*const.m]; %force net on body C2 frame basis
T = torque1(v0,w0,const.r0,const.h0,const.Lt,const.Ln); %torque net acting on body C2 frame basis
% 6dof system 13 equation variation
% cheat sheet for state-space
% x1 = x
% x2 = y
% x3 = z
% x4 = u
% x5 = v
% x6 = w
% x7 = wx
% x8 = wy
% x9 = wz
% x10 = e0
% x11 = e1
% x12 = e2
% x13 = e3

xdot =[
    x(4) - x(8)*x(3) + x(9)*x(2); %xdot
    x(5) + x(7)*x(3) - x(9)*x(1); %ydot
    x(6) - x(7)*x(2) + x(8)*x(1); %zdot
    f(1)/const.m - x(8)*x(6)  + x(9)*x(5); %udot
    f(2)/const.m + x(7)*x(6)  - x(9)*x(4); %vdot
    f(3)/const.m - x(7)*x(5)  + x(8)*x(4); %wdot
    T(1)/Ixx + (Iyy - Izz)*x(8)*x(9)/Ixx; %wxdot
    T(2)/Iyy + (Izz - Ixx)*x(7)*x(9)/Iyy; %wydot
    T(3)/Izz + (Ixx - Iyy)*x(7)*x(8)/Izz; %wzdot
    -2*x(11)*x(7)/4 - 2*x(12)*x(8)/4 - 2*x(13)*x(9)/4; %e0dot
    2*x(10)*x(7)/4 - 2*x(13)*x(8)/4 + 2*x(12)*x(9)/4; %e1dot
    2*x(13)*x(7)/4 + 2*x(10)*x(8)/4 - 2*x(11)*x(9)/4; %e2dot
    2*x(11)*x(8)/4 - 2*x(12)*x(7)/4 + 2*x(10)*x(9)/4; %e3dot
]
end
%%


















