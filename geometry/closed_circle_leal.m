%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y,high,low,alpha] = closed_circle_leal(x1,y1,N,DELTA)

theta = 0:2*pi/N:2*pi;
%second legendre polznomial

a = 1;
b = (1-DELTA)/(1+DELTA);

%build the circle (half)
X = a.*cos(theta)+x1;
Y = b.*sin(theta)+y1;

%rescale in order to have volume equal to 4/3*pi
V = axis_int(X,Y);
alpha = nthroot(4/3*pi/V,3);
a = alpha*a;
b = alpha*b;

%build the new circle (half)
X = a.*cos(theta)+x1;
Y = b.*sin(theta)+y1;

ds = sqrt((X(1:end-1)-X(2:end)).^2+(Y(1:end-1)-Y(2:end)).^2);
high=max(ds);
low=min(ds);

end