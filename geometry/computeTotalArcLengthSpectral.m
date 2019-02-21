%compute Arc Length Spectrally

function l = computeTotalArcLengthSpectral(x,y,SPECTRAL)

%drivative and integration
D1 = SPECTRAL.D1;
w = SPECTRAL.WG;

%deribvatives
xp = D1*x;  yp = D1*y;

%function to integrate
int = sqrt(xp.^2+yp.^2);

%perform integration
l = w'*int;