%compute Arc Length Spectrally

function l = computeArcLengthSpectral(x,y,SPECTRAL)

%drivative and integration
D1 = SPECTRAL.D1;

t = SPECTRAL.t;
xp = D1*x;  yp = D1*y;

int = sqrt(xp.^2+yp.^2);

dt = diff(t);
l = [0; cumsum((int(1:end-1)+int(2:end))/2.*dt)];