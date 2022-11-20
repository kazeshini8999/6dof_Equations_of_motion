function [f] = drag1(v0,w0,r0,h0,Lt,Ln,b)
%drag1 returns the drag force acting on the cylinder in terms of the state variables 
u = v0(1);
v = v0(2);
w = v0(3);
wx = w0(1);
wy = w0(2);
wz = w0(3);
%% 1 top plane
df1t_x = @(r,p) -Lt.*r.*((w - r.*wy.*cos(p) + r.*wx.*sin(p)).^2 + (v - (h0.*wx)./2 + r.*wz.*cos(p)).^2 + (u + (h0.*wy)./2 - r.*wz.*sin(p)).^2).^(b/2).*(u + (h0.*wy)./2 - r.*wz.*sin(p));
f1t_x = integral2(df1t_x,0,r0,0,2*pi,Method="auto");
df1t_y = @(r,p) -Lt.*r.*((w - r.*wy.*cos(p) + r.*wx.*sin(p)).^2 + (v - (h0.*wx)./2 + r.*wz.*cos(p)).^2 + (u + (h0.*wy)./2 - r.*wz.*sin(p)).^2).^(b./2).*(v - (h0.*wx)./2 + r.*wz.*cos(p));
f1t_y = integral2(df1t_y,0,r0,0,2*pi,Method="auto");
f1t_z = 0;
f1n_x = 0;
f1n_y = 0;
df1n_z = @(r,p) -Ln.*r.*sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*(sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)) + 1).*((w - r.*wy.*cos(p) + r.*wx.*sin(p)).^2 + (v - (h0.*wx)./2 + r.*wz.*cos(p)).^2 + (u + (h0.*wy)./2 - r.*wz.*sin(p)).^2).^(b./2).*(w - r.*wy.*cos(p) + r.*wx.*sin(p));
f1n_z = integral2(df1n_z,0,r0,0,2*pi,Method="auto");
f1 = [f1t_x + f1n_x;f1t_y + f1n_y;f1t_z+f1n_z]; %total force on top plane
%% 2 bottom plane
df2t_x = @(r,p) Lt.*r.*((w - r.*wy.*cos(p) + r.*wx.*sin(p)).^2 + (v + (h0.*wx)./2 + r.*wz.*cos(p)).^2 + ((h0.*wy)./2 - u + r.*wz.*sin(p)).^2).^(b./2).*((h0.*wy)./2 - u + r.*wz.*sin(p));
f2t_x  = integral2(df2t_x,0,r0,0,2*pi,Method="auto");
df2t_y = @(r,p) -Lt.*r.*(v + (h0.*wx)./2 + r.*wz.*cos(p)).*((w - r.*wy.*cos(p) + r.*wx.*sin(p)).^2 + (v + (h0.*wx)./2 + r.*wz.*cos(p)).^2 + ((h0.*wy)./2 - u + r.*wz.*sin(p)).^2).^(b./2);
f2t_y = integral2(df2t_y,0,r0,0,2*pi,Method="auto");
f2t_z = 0;
f2n_x = 0;
f2n_y = 0;
df2n_z = @(r,p) -Ln.*r.*sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*(sign(w - r.*wy.*cos(p) + r.*wx.*sin(p)) - 1).*(w - r.*wy.*cos(p) + r.*wx.*sin(p)).*((w - r.*wy.*cos(p) + r.*wx.*sin(p)).^2 + (v + (h0.*wx)./2 + r.*wz.*cos(p)).^2 + ((h0.*wy)./2 - u + r.*wz.*sin(p)).^2).^(b/2);
f2n_z = integral2(df2n_z,0,r0,0,2*pi,Method="auto");
f2 = [f2t_x + f2n_x;f2t_y + f2n_y;f2t_z+f2n_z]; %total force on bottom plane
%% 3 lateral force
df3n_x = @(h,p) -Ln.*r0.*sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)).*cos(p).*(sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)) + 1).*((w - r0.*wy.*cos(p) + r0.*wx.*sin(p)).^2 + (v - h.*wx + r0.*wz.*cos(p)).^2 + (u + h.*wy - r0.*wz.*sin(p)).^2).^(b./2).*(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p));
f3n_x = integral2(df3n_x,-h0/2,h0/2,0,2*pi,Method="auto");
df3n_y = @(h,p) -Ln.*r0.*sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)).*sin(p).*(sign(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p)) + 1).*((w - r0.*wy.*cos(p) + r0.*wx.*sin(p)).^2 + (v - h.*wx + r0.*wz.*cos(p)).^2 + (u + h.*wy - r0.*wz.*sin(p)).^2).^(b./2).*(u.*cos(p) + v.*sin(p) + h.*wy.*cos(p) - h.*wx.*sin(p));
f3n_y = integral2(df3n_y,-h0/2,h0/2,0,2*pi,Method="auto");
f3n_z = 0;
df3t_x = @(h,p) Lt.*r0.*((w - r0.*wy.*cos(p) + r0.*wx.*sin(p)).^2 + (v - h.*wx + r0.*wz.*cos(p)).^2 + (u + h.*wy - r0.*wz.*sin(p)).^2).^(b/2).*((u.*cos(2.*p))./2 - u./2 + (v.*sin(2.*p))./2 - (h.*wy)./2 + r0.*wz.*sin(p) + (h.*wy.*cos(2.*p))./2 - (h.*wx.*sin(2.*p))./2);
f3t_x = integral2(df3t_x,-h0/2,h0/2,0,2*pi,Method="auto");
df3t_y = @ (h,p) Lt.*r0.*cos(p).*((w - r0.*wy.*cos(p) + r0.*wx.*sin(p)).^2 + (v - h.*wx + r0.*wz.*cos(p)).^2 + (u + h.*wy - r0.*wz.*sin(p)).^2).^(b./2).*(u.*sin(p) - v.*cos(p) - r0.*wz + h.*wx.*cos(p) + h.*wy.*sin(p));
f3t_y = integral2(df3t_y,-h0/2,h0/2,0,2*pi,Method="auto",AbsTol=1e-09);
df3t_z = @ (h,p) Lt.*r0.*((w - r0.*wy.*cos(p) + r0.*wx.*sin(p)).^2 + (v - h.*wx + r0.*wz.*cos(p)).^2 + (u + h.*wy - r0.*wz.*sin(p)).^2).^(b./2).*((u.*cos(2.*p))./2 - u./2 + (v.*sin(2.*p))./2 - (h.*wy)./2 + r0.*wz.*sin(p) + (h.*wy.*cos(2.*p))./2 - (h.*wx.*sin(2.*p))./2);
f3t_z = integral2(df3t_z,-h0/2,h0/2,0,2*pi,Method="auto");
f3 = [f3t_x + f3n_x;f3t_y + f3n_y;f3t_z+f3n_z]; %total force on the lateral surface (needed only for 3d coin)
%%
f= f1+f2+f3; %force net return value
end

