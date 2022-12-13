function[x1] = impact_solver(t,x)
%function to solve newtons hypothesis alongside momentum conservation and
%constraint equations to return post impact state-space vector to feed into
%ode 89 again 
xi = transpose(x(length(t),1:end)); % pre impact state-space in C2 basis
C2cO = Euler2DCM(xi(10:13)); %rotation matrix to get to C2 space from O space
OcC2 = transpose(C2cO); % rotation matrix to get to O space from C2 space
psi = C2cO(3,3); % angle between Ko and Kc2 axes
phi = real(impact_point(t,x)); %angle to define the polar co-ords of the point of impact
w = OcC2*xi(7:9); %angular vel of C2 wrt to O in O frame basis pre impact
v = OcC2*xi(4:6); % vel of COM wrt to O in C2 space pre impact
rc2 = [const.r0*cos(phi(1));const.r0*sin(phi(1));-sign(psi)*const.h0/2];
r = OcC2*[const.r0*cos(phi(1));const.r0*sin(phi(1));-sign(psi)*const.h0/2]; %position vector of A wrt to C2 in C2 frame basis
va = real(OcC2*(vec2cross(xi(7:9))*rc2 + xi(4:6))); %vel of Point A wrt to O in O frame basis pre impact
B = [va(1);va(2);-va(3)*const.chi]; % vel of point A wrt to O in C2 frame basis post impact (restitution) 
I = OcC2*PInertiaTens(const.m,const.r0,const.h0)*C2cO; %inertia tensor of coin

%% defining the coeffs of the matrix A for the Ax = B sys of equations
a11 =1;
a18 = r(3);
a19 = -r(2);
b1 = B(1);

a22 = 1;
a29 = r(1);
a27 = -r(3);
b2 = B(2);

a33  = 1;
a37 = r(2);
a38 = -r(1);
b3 = B(3);

a44  = 0;
a55 = 0;
a66 = -1;
a41 = const.m;
a52 = const.m;
a63 = const.m;
b4 = const.m*v(1);
b5 = const.m*v(2);
b6 = const.m*v(3);

a77 = I(1,1);
a76 = -r(2);
a75 = r(3);

a84 = -r(3);
a86 = r(1);
a88 = I(2,2);

a94 = r(2);
a95 = -r(1);
a99 = I(3,3);

b7 = I(1,1)*w(1);
b8 = I(2,2)*w(2);
b9 = I(3,3)*w(3);



A = [a11 0 0 0 0 0 0 a18 a19;
    0 a22 0 0 0 0 a27 0 a29;
    0 0 a33 0 0 0 a37 a38 0;
    a41 0 0 a44 0 0 0 0 0;
    0 a52 0 0 a55 0 0 0 0;
    0 0 a63 0 0 a66 0 0 0;
    0 0 0 0 a75 a76 a77 0 0;
    0 0 0 a84 0 a86 0 a88 0;
    0 0 0 a94 a95 0 0 0 a99]; %main part, use symbolic toolbox in a different script to generate the coeffs of the matrix
b = transpose([b1 b2 b3 b4 b5 b6 b7 b8 b9]);%main part, use sym toolbox in a diff script to get values
xo = A\b; % moneyshot
velc2 = C2cO*[xo(1:3)];
wc2 = C2cO*[xo(7:9)];
x1 = real([xi(1);xi(2);xi(3);velc2(1);velc2(2);velc2(3);wc2(1);wc2(2);wc2(3);xi(10);xi(11);xi(12);xi(13)]);

end

%%
