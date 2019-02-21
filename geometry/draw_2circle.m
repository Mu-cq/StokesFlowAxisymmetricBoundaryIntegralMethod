%draw 2 ellipses

function [X,Y]=draw_2circle(x1,y1,N,DELTA1,DELTA2)

theta = 0:pi/N:pi;
%second legendre polznomial

a = 1;
b = (1-DELTA1)/(1+DELTA1);

%build the circle (half)
X = a.*cos(theta)+x1;
Y = b.*sin(theta)+y1;

%rescale in order to have volume equal to 4/3*pi
V = axis_int(X,Y);
alpha = nthroot(4/3*pi/V,3);
a = alpha*a;
b = alpha*b;

%build the new circle (half)
X1 = a.*cos(theta)+x1;
Y1 = b.*sin(theta)+y1;

a = 1;
b = (1-DELTA2)/(1+DELTA2);

%build the circle (half)
X = a.*cos(theta)+x1;
Y = b.*sin(theta)+y1;

%rescale in order to have volume equal to 4/3*pi
V = axis_int(X,Y);
alpha = nthroot(4/3*pi/V,3);
a = alpha*a;
b = alpha*b;

%build the new circle (half)
X2 = a.*cos(theta)+x1+2.5;
Y2 = b.*sin(theta)+y1;

X = [X1' X2'];
Y = [Y1' Y2'];

end