function[x1] = impact_solver(t,x)
%function to solve newtons hypothesis alongside momentum conservation and
%constraint equations to return post impact state-space vector to feed into
%ode 89 again 
xi = transpose(x(length(t),1:end)); % pre impact state-space in C2 basis
C2cO = Euler2DCM(xi(10:13)); %rotation matrix to get to C2 space from O space
OcC2 = transpose(C2cO); % rotation matrix to get to O space from C2 space
w_C2 = xi(7:9); %angular vel of C2 wrt to O in C2 frame basis
w = OcC2*w_C2; %angular vel of C2 wrt to O in O frame basis
v_C2 = xi(4:6); % vel of COM wrt to O in C2 space
v = OcC2*v_C2; % vel of COM wrt to O in O spcace
phi = impact_point(x(length(t),1:end)); %angle to define the polar co-ords of the point of impact
psi = OcC2(3,3); 
r_A_C2_C2 = [const.r0*cos(phi(1));const.r0*sin(phi(1));-sign(cos(psi))*const.h0/2];
r = OcC2*r_A_C2_C2;
Icm_c2 = PInertiaTens(const.m,const.r0,const.h0);
I = OcC2*Icm_c2*C2cO;
%% defining the coeffs of the matrix A for the Ax = B sys of equations
a11 =const.m;
a14 = -1;
b1 = const.m*v(1);
a22 = const.m;
a25 = -1;
b2 = const.m*v(2);
a33 = const.m;
a36 = -1;
b3 = const.m*v(3);
a41 = 1;
a52 = 1;
a63 = 1;
b4 = 0;
b5 = 0;
b6 = const.chi;
a77 = I(1,1);
a78 = I(1,2);
a79 = I(1,3);
a76 = -r(2);
a75 = r(3);
a84 = r(3);
a86 = -r(1);
a87 = I(1,2);
a88 = I(2,2);
a89 = I(2,3);
a94 = r(2);
a95 = -r(1);
a97 = I(1,3);
a98 = I(2,3);
a99 = I(3,3);

b7 = I(1,1)*w(1) + I(1,2)*w(2) + I(1,3)*w(3);
b8 = I(1,2)*w(1) + I(2,2)*w(2) + I(2,3)*w(3);
b9 = I(1,3)*w(1) + I(2,3)*w(2) + I(3,3)*w(3);



A = [a11 0 0 a14 0 0 0 0 0;
    0 a22 0 0 a25 0 0 0 0;
    0 0 a33 0 0 a36 0 0 0;
    a41 0 0 0 0 0 0 0 0;
    0 a52 0 0 0 0 0 0 0;
    0 0 a63 0 0 0 0 0 0;
    0 0 0 0 a75 a76 a77 a78 a79;
    0 0 0 a84 0 a86 a87 a88 a89;
    0 0 0 a94 a95 0 a97 a98 a99]; %main part, use symbolic toolbox in a different script to generate the coeffs of the matrix
b = transpose([b1 b2 b3 b4 b5 b6 b7 b8 b9]);%main part, use sym toolbox in a diff script to get values
xo = A\b; % moneyshot
xo1 = C2cO*[xo(1:3)];
xo2 = C2cO*[xo(7:9)];
%convert v and w into C2 space concat vector taking unchanged values and
%return function output
x1 = real([xi(1);xi(2);xi(3);xo1(1);xo1(2);xo1(3);xo2(1);xo2(2);xo2(3);xi(10);xi(11);xi(12);xi(13)]);


end

%%
