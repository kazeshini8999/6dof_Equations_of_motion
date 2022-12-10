hold on
% code snippet to plot center of mass trajectory
x3= zeros(size(t));
y3 = zeros(size(t));
z3 = zeros(size(t));
zdot = zeros(size(t));
for K = 1:size(t)
    rc1 = [x(K,1);x(K,2);x(K,3)];
    e = [x(K,10);x(K,11);x(K,12);x(K,13)];
    rc0 = transpose(Euler2DCM(e))*rc1;
    zdt = x(6) - x(K,7)*x(K,2) + x(K,8)*x(K,1);
    R = transpose(Euler2DCM(e));
    
    x3(K) = rc0(1);
    y3(K) = rc0(2);
    z3(K) = rc0(3);
    zdot(K) = x(K,6) - x(K,7)*x(K,2) + x(K,8)*x(K,1);
    ei(K) =  x(K,10)^2 + x(K,11)^2 + x(K,12)^2 + x(K,13)^2;
end
% plot(eii,t)
% xlim([0.995 1.005])
% e1= [x(450,10);x(450,11);x(450,12);x(450,13)];
% R = transpose(Euler2DCM(e1))
% quat2dcm(transpose(e1))

plot3(x3,y3,z3)
%plot3(vertcat(x1,x2),vertcat(y1,y2),vertcat(z1,z2))
%view(0,0)
% plot(t,zdot)!
grid on
xlabel('x')
ylabel('y')

zlabel('z')