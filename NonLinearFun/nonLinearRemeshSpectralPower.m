%non linear function for spectral remesh

function R = nonLinearRemeshSpectralPower(expon,t,x,y,SPECTRAL)

x = x(t.^expon);   y = y(t.^expon);

%compute arclenght
l = computeArcLengthSpectral(x,y,SPECTRAL);
%l = computeArcLength(x,y);
l0 = l(end);
l = l/l0;

%R = norm(l-t)/numel(x);
R = max(abs(l-t));