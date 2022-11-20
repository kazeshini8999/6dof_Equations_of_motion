function [e] = DCM2Euler(phi,theta,psi,i,j,k)
%DCM2Euler returns Euler parameters corresponding to a particular euler
%angle sequence
R = DCM(k,psi)*DCM(j,theta)*DCM(i,phi); %rotation matrix 
e0 = sqrt(trace(R)/4 + 1/4);
e1 = (R(3,2) - R(2,3))/(4*e0);
e2 = (R(1,3) - R(3,1))/(4*e0);
e3 = (R(2,1) - R(1,2))/(4*e0);
e = [e0;e1;e2;e3];
end