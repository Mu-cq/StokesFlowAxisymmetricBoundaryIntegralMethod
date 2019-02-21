%eqution for the circle

function [x,y] = circleEquation(theta,PARAM,panel)

R = PARAM.rArc(panel);
x0 = PARAM.x0_Circle(panel);
y0 = PARAM.y0_Circle(panel);

x = R*cos(theta) + x0;
y = R*sin(theta) + y0;