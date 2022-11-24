syms rx ry rz I11 I22 I33 I12 I13 I23 wx wy wz x7 x8 x9 x4 x5 x6
I = [I11 I12 I13;I12 I22 I23;I13 I23 I33];
eqn1 = vec2cross([rx;ry;rz])*[x4;x5;x6] - I*[(x7 - wx);(x8 - wy);(x9 - wz)];
expand(eqn1)

