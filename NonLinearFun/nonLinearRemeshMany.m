%non linear function for spectral remesh

function R = nonLinearRemeshMany(t,tOld,x,y,SPECTRAL)

if PARAM.legendre==1
    %t = [tOld(1); t(2:end-1); tOld(end)];
elseif SPECTRAL.legendre==0
    t = [tOld(1); t; tOld(end)];
end
x = x(t);   y = y(t);

%compute arclenght
l = computeArcLengthSpectral(x,y,SPECTRAL);
%l = computeArcLength(x,y);
l0 = l(end);
l = l/l0;

R = l-tOld;