%eqution for the circle

function [a,b,c] = lineCoeff(PARAM,panel)

xStart = PARAM.xStart(panel);
yStart = PARAM.yStart(panel);
xEnd = PARAM.xEnd(panel);
yEnd = PARAM.yEnd(panel);

m = (yEnd-yStart)/(xEnd-xStart);
c = yStart-m*xStart;
a = m;
b = -1;