%compute complex potential

clear all
close all

%parameters
Q1 = 0;
Q2 = [1 1];
%Q2 = ones(1,1000);

%phisical space
x1 = 0; y1 = 0;   %source_axis
x2 = 0; y2 = 0.1;   %source_in the domain
x = linspace(-5,5,100);
y = linspace(0,3,30);
[X0,Y0] = meshgrid(x,y);
THETA = atan(Y0./X0).*(X0>0) + (atan(Y0./X0)+pi).*(X0<0) + pi/2*(X0==0);

vx_ring = 0;
vy_ring = 0;

for i = 1:numel(Q2)
    X2 = repmat(x2,numel(y),numel(x));
    Y2 = repmat(y2,numel(y),numel(x));
    [SX,SY] = poisson_ax_fs (X2,Y2,X0,Y0);
    y2 = y2+2;
    
    %analytical formualtion for a ring point source
    vx_ring = vx_ring + Q2(i)/(4*pi)*SX;
    vy_ring = vy_ring + Q2(i)/(4*pi)*SY;
    
    %Q2 = Q2*y2;
end

%analytical formualtion for a point source in 3D
r = sqrt((X0-x1).^2+(Y0-y1).^2);
vr_point = Q1./(4*pi*r.^2);

%plot velocity field
figure
vx = vr_point.*cos(THETA)+vx_ring;  vy = vr_point.*sin(THETA)+vy_ring;
no = sqrt(vx.^2+vy.^2);
quiver(X0,Y0,vx./no,vy./no)
axis equal
title('velocity field')
%grid on

