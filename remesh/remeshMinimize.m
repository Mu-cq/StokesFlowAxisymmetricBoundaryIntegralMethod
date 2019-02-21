%non linear function for spectral remesh

function R = remeshMinimize(t,tOld,x,y,SPECTRAL)

D1 = SPECTRAL.D1;
w = SPECTRAL.WG';

if SPECTRAL.legendre==1
    %t = [tOld(1); t(2:end-1); tOld(end)];
elseif SPECTRAL.legendre==0
    t = [tOld(1); t(2:end-1); tOld(end)];
end
x = x(t);   y = y(t);

%compute arclenght
l = computeTotalArcLengthSpectral(x,y,SPECTRAL);

%compute metric term (first derivative)
dl = sqrt((D1*x).^2+(D1*y).^2)/l;
ddl = D1*dl;

%residuals
if SPECTRAL.normRemesh==1
    R(1) = w*(dl-(D1*SPECTRAL.remeshMapping(SPECTRAL.t))).^2;
elseif SPECTRAL.normRemesh==2
    error('Not updated for mapping')
    R(1) = w*ddl.^2;
elseif SPECTRAL.normRemesh==3
    R(1) = w*(dl-(D1*SPECTRAL.remeshMapping(SPECTRAL.t))).^2 + w*ddl.^2;
end