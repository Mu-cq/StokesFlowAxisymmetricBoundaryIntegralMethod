%non linear function for spectral remesh

function R = nonLinearRemeshSpectralPowerTwo(expon,t,x,y,SPECTRAL)

expon1 = expon(1);  expon2 = expon(2);
x = x(t.^expon1);   y = y(t.^expon2);

%compute arclenght
l = computeArcLengthSpectral(x,y,SPECTRAL);
%l = computeArcLength(x,y);
l0 = l(end);
l = l/l0;

R = norm(l-t)/numel(x);
%R = max(abs(l-t));