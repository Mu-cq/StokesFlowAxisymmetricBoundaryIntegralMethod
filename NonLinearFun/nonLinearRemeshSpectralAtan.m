%non linear function for spectral remesh

function R = nonLinearRemeshSpectralAtan(expon,t,x,y,SPECTRAL)

%x = x((t.^expon-0.5)-atan(-0.5));   y = y((t.^expon-0.5)-atan(-0.5));
x = (x(t-0.5)+atan(0.5))./2/atan(0.5)./expon;
y = (y(t-0.5)+atan(0.5))./2/atan(0.5)./expon;

%compute arclenght
l = computeArcLengthSpectral(x,y,SPECTRAL);
%l = computeArcLength(x,y);
l0 = l(end);
l = l/l0;

%R = norm(l-t)/numel(x);
R = max(abs(l-t));