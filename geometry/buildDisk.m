%build disk (rectangle in axisymmrtric)

function [x,y] = buildDisk(PARAM,res)

%get parameters
L = PARAM.L;
R = PARAM.R;
x0 = PARAM.x0;

%number of point
m = L*res;  n = res*R;

%build disk
x1 = x0-L/2;    y1 = 0;
x2 = x0-L/2;    y2 = R;
[xLeft,yLeft] = drawline2(x1,y1,x2,y2,n);

x1 = x0-L/2;    y1 = R;
x2 = x0+L/2;    y2 = R;
[xUp,yUp] = drawline2(x1,y1,x2,y2,m);

x1 = x0+L/2;    y1 = R;
x2 = x0+L/2;    y2 = 0;
[xRight,yRight] = drawline2(x1,y1,x2,y2,n);

x = [xLeft(1:end-1) xUp(1:end-1) xRight];
y = [yLeft(1:end-1) yUp(1:end-1) yRight];