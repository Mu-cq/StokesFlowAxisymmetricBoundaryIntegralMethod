%compute Arc Length Spectrally

function l = computeArcLengthSpectral2(xy,SPECTRAL)

x = xy(1:2:end-1);
y = xy(2:2:end);

%drivative and integration
D1 = SPECTRAL.D1;

t = SPECTRAL.t;
xp = D1*x;  yp = D1*y;

%integration
%l = zeros(numel(x),1);
int = sqrt(xp.^2+yp.^2);
%l0 = w'*[int(1); zeros(numel(x)-1,1)];
% for i = 1:numel(x)
%     
%     l(i) = w'*[int(1:i); zeros(numel(x)-i,1)]-l0;
%     
% end

%l = [0; cumsum(sqrt(diff(x).^2+diff(y).^2))];
dt = diff(t);
l = [0; cumsum((int(1:end-1)+int(2:end))/2.*dt)];