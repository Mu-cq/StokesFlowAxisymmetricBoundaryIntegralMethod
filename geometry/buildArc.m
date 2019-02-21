%build arc

function [x,y] = buildArc(r,x0,y0,inOut,PARAM,nElem)

n = PARAM.n(nElem)+1;

theta = linspace(inOut(1),inOut(end),n);
x = r*cos(theta) + x0;
y = r*sin(theta) + y0;