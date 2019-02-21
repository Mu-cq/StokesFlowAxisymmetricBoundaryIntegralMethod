%compute channel wall with potential theory in 2D, wall is the x axis

clear all
close all

%physical variables
y0 = -1;     %source position
x = -5:0.2:5;
y  = -1:0.1:1;
[X,Y] = meshgrid(x,y);
Q = 10;  %flow rate of injection

%from analytical form of the complex potential
den = (X.^2-Y.^2+y0^2).^2+4*X.^2.*Y.^2;
num_x = X.^3+X*y0^2+X.*Y.^2;
num_y = -Y.^3+y0^2*Y-2*X.^2.*Y;
%VELOCITIES
u = Q/pi*num_x./den;    %velocity in x direction
v = Q/pi*num_y./den;     %velocity in y direction
%VELOCITIES FIRST DERIVATIVES
dudx = ((3*X.^2+Y.^2+y0^2).*den-num_x.*(4*X.^3+4*X.*Y.^2+4*X*y0^2))./den.^2;


figure
quiver(X,Y,u,v)
axis equal
title('velocity field')
xlabel('x')
ylabel('y')

figure
contour(X,Y,dudx,500)
axis equal
title('\partial u / \partial x')
xlabel('x')
ylabel('y')