%find inverse of a 1D function

function myNewXY = findInv(x,y,t,SPECTRAL)

%initial guess
xy = zeros(2*numel(x),1);
xy(1:2:end-1) = x;
xy(2:2:end) = y;

%find inverse
f = @(xyUnk) max(abs(computeArcLengthSpectral2(xyUnk,SPECTRAL)-t));
%options = optimoptions('fsolve','TolFun',1e-10,'TolX',1e-10,'Display','iter','Algorithm','Levenberg-Marquardt');
options = optimoptions('fsolve','TolFun',1e-10,'TolX',1e-10,'Display','iter');
myNewXY = fsolve(f,xy,options);