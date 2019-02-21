%compute radial vector with spectral methods

function [nx,ny] = radialVectorSpectral(x,y,SPECTRAL)

%compute center of mass
xcm = CenterMassCurvAxisSpectral(x,y,SPECTRAL);

%vector norm
norm = sqrt((x-xcm).^2+y.^2);

nx = (x-xcm)./norm;
ny = y./norm;