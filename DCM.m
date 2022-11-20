function [R] = DCM(i,l)
% takes in one integer from [0,3] defining which axis the rotation is
% happening around, and one symbolic/numeric variable defining the angle
% rotated.
    function d = kronDel(j,k)

        if j == k
            d = 1;
        else
            d = 0;
        end
    end
    r1 = [1 0 0;0 cos(l) sin(l);0 -sin(l) cos(l)];
    r2 = [cos(l) 0 -sin(l);0 1 0;sin(l) 0 cos(l)];
    r3 = [cos(l) sin(l) 0;-sin(l) cos(l) 0;0 0 1];
    R = kronDel(i,1)*r1 + kronDel(i,2)*r2 + kronDel(i,3)*r3;
    
end

