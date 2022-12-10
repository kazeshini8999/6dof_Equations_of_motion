function [e] = DCM2Euler(phi,theta,psi)
%DCM2Euler returns Euler parameters corresponding to a particular euler
%angle sequence
R = DCM(3,psi)*DCM(2,theta)*DCM(1,phi); %rotation matrix 
e0 = sqrt(trace(R)/4 + 1/4);
e1 = (R(3,2) - R(2,3))/(4*e0);
e2 = (R(1,3) - R(3,1))/(4*e0);
e3 = (R(2,1) - R(1,2))/(4*e0);
e = [e0;e1;e2;e3];
end