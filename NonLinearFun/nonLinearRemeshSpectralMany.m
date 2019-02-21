%non linear function for spectral remesh

function R = nonLinearRemeshSpectralMany(t,tOld,x,y,SPECTRAL)

D1 = SPECTRAL.D1;

t = [tOld(1); t; tOld(end)];
x = x(t);   y = y(t);

%compute arclenght
l = computeArcLengthSpectral(x,y,SPECTRAL);
%l = computeArcLength(x,y);
l0 = l(end);
l = l/l0;

R = D1*l-1;