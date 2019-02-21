%non linear function for spectral remesh

function R = nonLinearRemeshSpectralMany2(t,tOld,x,y,SPECTRAL)

D1 = SPECTRAL.D1;
%D2 = SPECTRAL.D2;

if SPECTRAL.legendre==1
    t = [tOld(1); t(2:end-1); tOld(end)];
elseif SPECTRAL.legendre==0
    t = [tOld(1); t; tOld(end)];
end
x = x(t);   y = y(t);

%compute arclenght
l = computeTotalArcLengthSpectral(x,y,SPECTRAL);

%compute metric term (first derivative)
dl = sqrt((D1*x).^2+(D1*y).^2)/l;

%residuals
R = dl-1;