%resize shape such that irt respect force free condition

function [xNew,yNew] = resizeShapeForceFree2(x,y,gamma)

%compute drop volume
V = axis_int_gauss_vect(x,y);

% from cartesian to polar coordinates
r = sqrt(x.^2+y.^2);
theta = atan(y./x);
theta = theta + pi*(theta<0);
    
%compute resizing factor
fStressX = @(unk) rescaleFunctionShape2(theta,r,gamma,unk,0,V);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','iter');
coeffs = fsolve(fStressX,[1 .01],options);

C0 = coeffs(1);
C2 = coeffs(2);

%legendre function
PPP = legendre(0,cos(theta));
P0 = PPP(1,:)';   %legendre polynomia
PPP = legendre(40,cos(theta));
P2 = PPP(1,:)';   %legendre polynomia

%new radius
r = r + C0*P0' + C2*P2';

%compute new droplet corrdinates
xNew = r.*cos(theta);
yNew = r.*sin(theta);