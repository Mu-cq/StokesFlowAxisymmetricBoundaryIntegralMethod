%non linear function for spectral remesh

function R = remeshSolve(t,tOld,x,y,l0,SPECTRAL)

D1 = SPECTRAL.D1;

if SPECTRAL.legendre==0||SPECTRAL.legendre==2
    t = [tOld(1); t(2:end-1); tOld(end)];
end

x = x(t);   y = y(t);

%compute arclenght
l = computeTotalArcLengthSpectral(x,y,SPECTRAL);

%compute metric term (first derivative)
xp = D1*x;  yp = D1*y;
dl = sqrt(xp.^2+yp.^2)/l;
ddl = D1*dl;

%residuals
if SPECTRAL.normRemesh==1
    %R = dl-objective;
    R = dl-(D1*SPECTRAL.remeshMapping(SPECTRAL.t));
    if SPECTRAL.legendre==1
        R = [R; l-l0];
    end
elseif SPECTRAL.normRemesh==2
    error('Not updated with mapping')
    R = [ddl; l-l0];
elseif SPECTRAL.normRemesh==3
    R = [(abs(dl-(D1*SPECTRAL.remeshMapping(SPECTRAL.t))) + abs(ddl)); l-l0];
end