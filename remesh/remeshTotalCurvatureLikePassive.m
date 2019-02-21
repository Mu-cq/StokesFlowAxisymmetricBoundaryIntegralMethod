%remesh clustering points in the regions oh huigh total curvature

function [x,y] = remeshTotalCurvatureLikePassive(x,y,howMuch)

disp('Remesh based on curvature')

%compute total curvature
[K1,K2] = computeCurvatureSplines(x',y',1);
K = K1.^2 + K2.^2 +0.005;

%compute arc lenght
dx = diff(x');   dy = diff(y');
dl = sqrt(dx.^2+dy.^2);
l = [0 cumsum(dl)];

%compute new elements size
dl = (1./K).^howMuch;
dl = (dl(1:end-1)+dl(2:end))/2;
lNew = [0 cumsum(dl)];
lNew = lNew/lNew(end)*l(end);

%new points
x = spline(l,x',lNew)';
y = spline(l,y',lNew)';