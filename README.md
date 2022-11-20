# 6dof_Equations_of_motion
6dof equations of motions to simulate the motion of a flipped coin using MatLab

Standard 13 equations 13 unknown system of equations are used, with the orientation of the coin described by euler-parameters instead of euler angles
The aerodynamic drag is modeled as a linear function of velocity (stokes law), however quadratic velocity dependance can also be used.
The solver used is ODE89. Simple Restitution with perfectly rough surface is used to model the collision. Matlab Poseplot is used to visualize the motion of the coin
