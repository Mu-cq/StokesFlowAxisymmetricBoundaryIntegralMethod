%pload initial shape form continuation method

function [x,y,xBase,yBase,perturb,PARAMminimalSeed] = uploadFromMinimalSeed(PARAM)

%load data
here = pwd;
cd(PARAM.res)
fileLoad = ['minimalSeedSpectralXYmodes_n=' num2str(PARAM.elemUpload) '_Ca=' num2str(PARAM.CaUpload) '_visc=' num2str(PARAM.viscUpload) '_Legendre=' num2str(PARAM.legendreUpload) '_BC=' num2str(PARAM.BCupload) '_T=' num2str(PARAM.ThorizonUpload) '_A0=' num2str(PARAM.A0upload) '.mat'];
minimalResults = load(fileLoad);
cd(here)
PARAMminimalSeed = minimalResults.PARAM;

%modes and perturbation of the base shape, this is useful for stability
%analysis
xBase = minimalResults.xBase;
yBase = minimalResults.yBase;
perturb = minimalResults.perturb;

%stability works well only with the same spectral grid and same number of
%harmonics
if PARAM.algorithm==5 && (PARAM.legendre~=PARAMminimalSeed.legendre)
    error('Spectral basis must be the same as for the computation of the base state')
elseif PARAM.algorithm==5 && (PARAM.dealiasing~=PARAMminimalSeed.dealiasing)
    error('Number of harmonics must be the same as in the base shape')
end

%compute modes from grid
x = minimalResults.x;
y = minimalResults.y;
[xMode,yMode] = fromGridToModes(x,y,PARAMminimalSeed);

%compute grid points on the new grid
if PARAMminimalSeed.legendre==1||PARAMminimalSeed.legendre==2
    
    x = interpLegendreZeroOne(PARAM.t,xMode);
    y = interpLegendreZeroOne(PARAM.t,yMode);
    
elseif PARAMminimalSeed.legendre==0
    
    x = chebcoeffs2chebvals(xMode);
    y = chebcoeffs2chebvals(yMode);
    
    x = chebfun(x,[0 1]);
    y = chebfun(y,[0 1]);
    
    x = x(PARAM.t);
    y = y(PARAM.t);
    
end












