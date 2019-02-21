%plot solid of revolution

function [ux,uy,uz] = velocityRevolution(X,Y,z,r,ur,uz)

theta = atan(Y./X) + pi*(X<0&&Y<0) + pi/2*(X<0&&Y>0);
R = sqrt(X.^2+Y.^2);

nGrid = 100;

%theta start end end
t1 = theta(1);
t2 = theta(2);

theta = linspace(t1,t2,nGrid);
theta = repmat(theta',1,numel(z));
z = repmat(z,nGrid,1);
r = repmat(r,nGrid,1);

%rotate shape
x = r.*cos(theta);
y = r.*sin(theta);

%plot3(x,y,z,'k')
surf(z,y,x)