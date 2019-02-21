%non linear function for spectral remesh

function R = nonLinearRemeshSpectral(tAdd,i,t,x,y,SPECTRAL)

tOld = t(i);
t(i) = tAdd;
x = x(t);   y = y(t);

%compute arclenght
l = computeArcLengthSpectral(x,y,SPECTRAL);
%l = computeArcLength(x,y);
l0 = l(end);
l = l/l0;

R = l(i)-tOld;