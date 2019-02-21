% constraint function for arc lenght

function [Disequality,Equality,jacobianDisequality,jacobianEquality] = constraintArcLenght(t,x,y,l0,SPECTRAL)

%current coordinate
x = x(t);   y = y(t);

%compute arclenght
l = computeTotalArcLengthSpectral(x,y,SPECTRAL);

D1 = SPECTRAL.D1;   D2 = SPECTRAL.D2;   w = SPECTRAL.WG';
xp = D1*x;  yp = D1*y;  xpp = D2*x;  ypp = D2*y;
h = sqrt(xp.^2+yp.^2);

J = w'.*((xp.*xpp+yp.*ypp)./h);

%residuals
Disequality = [];
jacobianDisequality = [];
Equality = l-l0;
jacobianEquality = h;