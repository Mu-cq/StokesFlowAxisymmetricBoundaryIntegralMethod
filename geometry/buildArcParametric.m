%build arc

function [x,y] = buildArcParametric(r,x0,y0,theta)

x = r*cos(theta) + x0;
y = r*sin(theta) + y0;