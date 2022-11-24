function [v] = cross2vec(m)
%cross2vec takes in a crossproduct matrix and returns the vector that
%formed it
v = [m(3,2);m(1,3);m(2,1)];
end