%build arc

function [x,y] = buildEllipse(x0,y0,inOut,PARAM,nElem,D)

%minor and major axis
b = ((1.0-D)/(1.0+D))^(1/3);
a = 1.0/b^2;

n = PARAM.n(nElem)+1;

theta = linspace(inOut(1),inOut(end),n);
x = a*cos(theta) + x0;
y = b*sin(theta) + y0;