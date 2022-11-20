function [m] = vector2cross(v)
%vec2cross takes a vector and returns the matrix that emulates a cross
%product in the same frame basis
m = [0 -v(3,1) v(2,1);v(3,1) 0 -v(1,1);-v(2,1) v(1,1) 0];
end