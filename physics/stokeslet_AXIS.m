%compute Stokeslet in axisymmetric domain, a velocity field in the Stokes flow due to a point
%force

clear all
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%point force location and direction
x0 = 0; y0 = 1.01; theta = 11*pi/3;

%create grid
x = -1:0.2:1;   y = 0:0.2:2;
[X,Y] = meshgrid(x,y);

%for vectorial operations
X0 = repmat(x0,numel(y),numel(x));  Y0 = repmat(y0,numel(y),numel(x));  THETA = repmat(theta,numel(y),numel(x));

%compute green's function
[SXX,SXY,SYX,SYY] = sgf_ax_fs_vect3 (X0,Y0,X,Y);

%compute velocity field
U = SXX.*cos(THETA) + SXY.*sin(THETA);
V = SYX.*cos(THETA) + SYY.*sin(THETA);

figure
quiver(X,Y,U,V)
hold on
quiver(x0,y0,cos(theta),sin(theta),'r')
plot(x0,y0,'or')
hold off
xlabel('x')
ylabel('y')