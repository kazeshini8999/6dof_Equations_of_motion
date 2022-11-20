function [s] = mag(v)
%mag returns the magnitude of a vector without the abs signs 
%   Detailed explanation goes here
s = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
end