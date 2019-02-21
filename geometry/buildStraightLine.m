%build arc

function [x,y] = buildStraightLine(x1,y1,x2,y2,PARAM,nElem)

n = PARAM.n(nElem)+1;

x = linspace(x1,x2,n);
y = linspace(y1,y2,n);