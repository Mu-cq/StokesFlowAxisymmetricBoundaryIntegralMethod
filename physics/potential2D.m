%plot potential flow field due to singular mass injection

%physical variables
Q = 1;
x0 = 0;  y0 = 0;    %injection point
x = linspace(-1,1,10);  y = linspace(-1,1,10);
[X,Y] = meshgrid(x,y);  %field points

%velocities and stresses computation
[u,v,nxx,nxy,nyx,nyy] = gf2D_laplace_fs(X,Y,x0,y0,Q);

%plot velocity field
figure
quiver(X,Y,u,v)
axis equal
xlabel('x')
ylabel('y')

%plot stresses
surf(X,Y,nyy)