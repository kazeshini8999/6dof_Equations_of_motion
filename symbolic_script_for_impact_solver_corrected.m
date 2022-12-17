syms x7 x8 x9 rx ry rz x1 x2 x3 x4 x5 x6 b1 b2 b3 u  v w wx wy wz m Ixx Iyy Izz Ixy Ixz Iyz
I = [Ixx Ixy Ixz;Ixy Iyy Iyz; Ixz Iyz Izz]
eqn1 = vec2cross([x6;x7;x8])*[rx;ry;rz] + [x3;x4;x5]-[x1;x2;b1]
eqn2 = [0;0;x9] -m*[x3;x4;x5] + m*[u;v;w]
eqn3 = vec2cross([rx;ry;rz])*[0;0;x9] - I*[x6;x7;x8] + I*[wx;wy;wz]