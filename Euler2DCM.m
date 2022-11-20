function [R] = Euler2DCM(e)
%Euler2DCM takes in Euler parameter/normalized quaternions and returns the
%rotation matrix to get to the inertial frame basis from body frame basis
e0 = e(1);
e1 = e(2);
e2 = e(3);
e3 = e(4);
R = [(e0^2 + e1^2 -e2^2-e3^2) (2*(e1*e2 + e0*e3)) (2*(e1*e3-e0*e2));
    (2*(e1*e2 - e0*e3)) (e0^2 -e1^2 + e2^2 - e3^2) (2*(e2*e3 + e0*e1));
    (2*(e1*e3 + e0*e2)) (2*(e2*e3 -e0*e1)) (e0^2-e1^2-e2^2+e3^2)];
end




