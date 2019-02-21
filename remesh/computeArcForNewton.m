%compute arc lenght spacing

function R = computeArcForNewton(xy,lMesh)

x = xy(1:2:end-1);
y = xy(2:2:end);

dx = diff(x);   dy = diff(y);   dl = sqrt(dx.^2+dy.^2);
l = [0; cumsum(dl)];

R = l-lMesh;