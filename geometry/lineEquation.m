%eqution for the circle

function [x,y] = lineEquation(t,PARAM,panel)

xStart = PARAM.xStart(panel);
yStart = PARAM.yStart(panel);
xEnd = PARAM.xEnd(panel);
yEnd = PARAM.yEnd(panel);

x = xStart*(1-t) + xEnd*t;
y = yStart*(1-t) + yEnd*t;