%rescale a shape to get a desired integral value

function R = rescaleFunctionShape2(theta,r,gamma,alpha,intFinal,V0)

%modify radius
%legendre function
PPP = legendre(0,cos(theta));
P0 = PPP(1,:)';   %legendre polynomia
PPP = legendre(40,cos(theta));
P2 = PPP(1,:)';   %legendre polynomia

%new radius
r = r + alpha(1)*P0' + alpha(2)*P2';
x = r.*cos(theta);  y = r.*sin(theta);

%compute stresses in the axial direction
dfX = normalStresses(x,y,gamma);

%compute residuals
R(1) = int_axis_spline_symmetric(x,y,dfX')-intFinal;
R(2) = axis_int_gauss_vect(x,y)-V0;