%plot species concentration for a Janus Particle

function [c,vSlip] = plotSpeciesFieldAnalyticalFunSmooth2(r,theta)

%k
k0 = 0.5;
k1 = 0.5;

%finite diff
D1 = finiteDifference1D(numel(theta),[2 0],1);

%coeff
cp0 = -k0/r;
cp1 = -k1/2/r^2;

%legendre poly
P1 = legendre(1,cos(theta));
P1 = P1(1,:)';
%P1der = D1*P1/pi;
P1der = -ones(numel(theta),1).*sin(theta);

%concentration
c = cp0 + cp1*P1;

%slip velocity
vSlip = cp1*P1der/r;