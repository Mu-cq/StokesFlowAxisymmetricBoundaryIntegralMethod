%compute Arc Length Spectrally

function l = computeArcLength(x,y)

dx2 = diff(x).^2;
dy2 = diff(y).^2;
int = sqrt(dx2+dy2);

l = [0; cumsum(int)];