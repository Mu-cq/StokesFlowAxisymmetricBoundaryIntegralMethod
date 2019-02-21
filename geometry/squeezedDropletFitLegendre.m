%create droplet shape in confined cone having continuous curvature

function [xDrop,yDrop,x,y] = squeezedDropletFitLegendre(R,alpha,h,thickness,theta,x0,L,nDrop)

%number of modes
modes = 15;

%parameters
PARAM.alpha = alpha;
PARAM.h = h;
PARAM.thickness = thickness;
PARAM.theta = theta;
PARAM.L = L;
PARAM.R = R;
PARAM.start = x0;
elem = 400;

%build the bubble interface at the first iteration
PARAM.Vtot = 4/3*pi*R^3*PARAM.alpha^3;
PARAM.h = PARAM.h + PARAM.thickness;
[R1,R2,H] = findR1R2H(PARAM);
[x,y] = draw_test3(0,0,H,R1,R2,elem,PARAM.theta);

%polar coordinates
r = sqrt(x.^2+y.^2);
theta = atan(y./x);
theta = theta + pi*(theta<0);

thetaGrid = linspace(0,pi,nDrop+1);

%fit legendre
f = LegendreSeriePolar(theta',r',modes,0);
rBuilt = LegendreBuildShape(thetaGrid,f,0);

%resize
errVhandle = @(dr) (axis_int_gauss_vect((rBuilt+dr)'.*cos(thetaGrid),(rBuilt+dr)'.*sin(thetaGrid))-PARAM.Vtot)/PARAM.Vtot;
dr = myNewtonMethodV2(errVhandle,0,1e-14,0);
rBuiltNew = rBuilt+dr;

%x,y coordinates
xDrop = rBuiltNew'.*cos(thetaGrid);
yDrop = rBuiltNew'.*sin(thetaGrid);

% figure
% plot(x,y)
% hold on
% plot(xDrop,yDrop)
% axis equal
% grid on









