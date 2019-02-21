%choose non linear function, 1 drop, Spectral

function [initial,x,y,Abase] = initialConditionDropMinimalSeed(perturb,V0,PARAM)

%upload shape
display('Upload initial shape from single newton method')
if PARAM.uploadShape==2
    [xGrid,yGrid] = uploadFromNewton(PARAM);
else
    error('You need to start from a base state obtained with Newton method')
end

%compute surface area of base state
Abase = surfaceCurvilinearAxisSpectral(xGrid,yGrid,PARAM);

%volume-conserving shape
[x,y] = getVolumeConservingShapeSpectral(perturb,xGrid,yGrid,V0,PARAM);

%place shape in the center of mass
if PARAM.placeShapeXCM==1
    
    display('Place shape center of mass in the origin')
    xcm = CenterMassCurvAxisSpectral(x,y,PARAM);
    x = x-xcm;
    
end

%go from grid points to modes
[xMode,yMode] = fromGridToModes(x,y,PARAM);
xMode(PARAM.dealiasing+1:end) = 0;  yMode(PARAM.dealiasing+1:end) = 0;
initial = zeros(2*PARAM.n+2,1);
initial(1:2:end-1) = xMode;     initial(2:2:end) = yMode;

if PARAM.remeshStart==1
    %execute remesh
    display('Remesh spectrally')
    initial = remeshDropSpectral(initial,PARAM);
end

