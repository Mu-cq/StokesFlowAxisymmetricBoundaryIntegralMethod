%test rigid body rotation

function [x1,y1] = rotatePoint(x0,y0,theta)

%matrix of rigid rotation
matr = [cos(theta) -sin(theta); sin(theta) cos(theta)];

%rotation
in = [x0; y0];
out = matr*in;
x1 = out(1);
y1 = out(2);









