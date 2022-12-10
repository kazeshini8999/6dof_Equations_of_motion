%code snippet to visualize the coin dynamics (rotation)

hold on
grid on

j = 681;
i = 200;
rc1 = [x(j,1);x(j,2);x(j,3)];
e = [x(j,10);x(j,11);x(j,12);x(j,13)];
q = quaternion(e(1),e(2),e(3),e(4));
rc0 = transpose(Euler2DCM(e))*rc1;
position  = [x3(j) y3(j) z3(j)]; % note the z component is scaled down by a factor of 100 to fit the coin travel and coin rotation in one plot
poseplot(q,position,"ENU",MeshFileName="cylinder1.stl",ScaleFactor=100);

% rc2 = [x(i,1);x(i,2);x(i,3)];
% e2 = [x(i,10);x(i,11);x(i,12);x(i,13)];
% q2 = quaternion(e2(1),e2(2),e2(3),e2(4));
% rc3 = transpose(Euler2DCM(e))*rc2;
% position2  = [x2(i) y2(i) z2(i)];
% poseplot(q2,position2,"ENU",MeshFileName="cylinder1.stl",ScaleFactor=100);
% view(90,0)

axis('equal')
hold off

    






%x1 = cos(u)*sin(v);
% y2 = sin(u)*sin(v);
% z3 = cos(v)*sin(v);
% % fsurf(x1,y2,z3)
% % axis equal
% % xyznew =R*[x1;y2;z3];
% % fsurf(xyznew(1), xyznew(2), xyznew(3))
% % xlabel('x')
% % ylabel('y')