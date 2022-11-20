function [t] = torque1(v0,w0,r0,h0,Lt,Ln)
%torque1 returns the torque acting on the cylinder as a function of state
%variables 
u = v0(1);
v = v0(2);
w = v0(3);
wx = w0(1);
wy = w0(2);
wz = w0(3);
%% 1 %all integrands are defined as function handles in vars r p and h
dT1_x = @(p,r) (Lt.*h0.*r.*(v - (h0.*wx)./2 + r.*wz.*cos(p)))./2 - Ln.*r.^2.*sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*sin(p).*(sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)) + 1).*(w - r.*wy.*cos(p) + r.*wx.*sin(p));
dT1_y = @(p,r) Ln.*r.^2.*sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*cos(p).*(sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)) + 1).*(w - r.*wy.*cos(p) + r.*wx.*sin(p)) - (Lt.*h0.*r.*(u + (h0.*wy)./2 - r.*wz.*sin(p)))./2;
dT1_z = @(p,r) Lt.*r.^2.*sin(p).*(u + (h0.*wy)./2 - r.*wz.*sin(p)) - Lt.*r.^2.*cos(p).*(v - (h0.*wx)./2 + r.*wz.*cos(p));
%% 2
dT2_x = @(p,r) -(Lt.*h0.*r.*(v + (h0.*wx)./2 + r.*wz.*cos(p)))./2 - Ln.*r.^2.*sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*sin(p).*(sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)) - 1).*(w - r.*wy.*cos(p) + r.*wx.*sin(p));
dT2_y = @(p,r) Ln.*r.^2.*sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*cos(p).*(sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)) - 1).*(w - r.*wy.*cos(p) + r.*wx.*sin(p)) - (Lt.*h0.*r.*((h0.*wy)./2 - u + r.*wz.*sin(p)))./2;
dT2_z = @(p,r) -Lt.*r.^2.*cos(p).*(v + (h0.*wx)./2 + r.*wz.*cos(p)) - Lt.*r.^2.*sin(p).*((h0.*wy)./2 - u + r.*wz.*sin(p));
%% 3
dT3_x = @(h,p) -h.*(Lt.*r0.*cos(p).*(u.*sin(p) - v.*cos(p) - r0.*wz + h.*wx.*cos(p) + h.*wy.*sin(p)) - Ln.*r0.*sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)).*sin(p).*(sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)) + 1).*(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p))) - Lt.*r0.^2.*sin(p).*(w - r0.*wy.*cos(p) + r0.*wx.*sin(p)); % x component of dtorque
dT3_y = @(h,p) h.*(Lt.*r0.*((u.*cos(2.*p))./2 - u./2 + (v.*sin(2.*p))./2 - (h.*wy)./2 + r0.*wz.*sin(p) + (h.*wy.*cos(2.*p))./2 - (h.*wx.*sin(2.*p))./2) - Ln.*r0.*sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)).*cos(p).*(sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)) + 1).*(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p))) + Lt.*r0.^2.*cos(p).*(w - r0.*wy.*cos(p) + r0.*wx.*sin(p));
dT3_z = @(h,p) Lt.*r0.^2.*(u.*sin(p) - v.*cos(p) - r0.*wz + h.*wx.*cos(p) + h.*wy.*sin(p));
%% integration to get total torques
T1_x = integral2(dT1_x,0,2*pi,0,r0,Method="auto",AbsTol=1e-05);
T1_y = integral2(dT1_y,0,2*pi,0,r0,Method="auto",AbsTol=1e-05);
T1_z = integral2(dT1_z,0,2*pi,0,r0,Method="auto",AbsTol=1e-05);
T1 = [T1_x;T1_y;T1_z];
%%
T2_x = integral2(dT2_x,0,2*pi,0,r0,Method="auto",AbsTol=1e-05);
T2_y = integral2(dT2_y,0,2*pi,0,r0,Method="auto",AbsTol=1e-06);
T2_z = integral2(dT2_z,0,2*pi,0,r0,Method="auto",AbsTol=1e-05);
T2 = [T2_x;T2_y;T2_z];
%%
T3_x = integral2(dT3_x,-h0/2,h0/2,0,2*pi,Method="auto",AbsTol=1e-05);
T3_y = integral2(dT3_y,-h0/2,h0/2,0,2*pi,Method="auto",AbsTol=1e-05);
T3_z = integral2(dT3_z,-h0/2,h0/2,0,2*pi,Method="auto",AbsTol=1e-05);
T3 = [T3_x;T3_y;T3_z];
%% return value
t = T1+T2+T3;
end