%pload initial shape form continuation method

function [x,y,newCa,xBase,yBase,perturb,PARAMnewton] = uploadFromContinuation(PARAM)

%load data
here = pwd;
cd(PARAM.res)
fileLoad = ['newtonMethodSpectralXYmodesCont_n=' num2str(PARAM.elemUpload) '_CaUp=' num2str(PARAM.CaUpUpload) '_CaDown='  num2str(PARAM.CaDownUpload) '_visc=' num2str(PARAM.viscUpload) '_Legendre=' num2str(PARAM.legendreUpload) '_BC=' num2str(PARAM.BCupload) '.mat'];
newtonResults = load(fileLoad);
cd(here)
PARAMnewton = newtonResults.PARAM;

%find derired shape
manyD = max(newtonResults.manyA)';
manyCa = newtonResults.manyCa;
xxx = newtonResults.manyA;
yyy = newtonResults.manyB;
[minVal,indMin] = min(sqrt(abs((manyD-PARAM.Dupload)/PARAM.Dupload).^2+abs((manyCa-PARAM.CaUpload)/PARAM.CaUpload).^2));
newCa = manyCa(indMin);

if minVal>1e-1
    warning('Upload shape is far from the desired one')
end

%grid points
x = xxx(:,indMin);
y = yyy(:,indMin);

if PARAM.algorithm==5
    error('Not implemented, it is better to upload from newton method if you want to perform stability analysis')
end

%modes and perturbation of the base shape, this is useful for stability
%analysis
xBase = newtonResults.xBase;
yBase = newtonResults.yBase;
perturb = newtonResults.perturb;

%stability works well only with the same spectral grid and same number of
%harmonics
if PARAM.algorithm==5 && (PARAM.legendre~=PARAMnewton.legendre)
    error('Spectral basis must be the same as for the computation of the base state')
elseif PARAM.algorithm==5 && (PARAM.dealiasing~=PARAMnewton.dealiasing)
    error('Number of harmonics must be the same as in the base shape')
end

%compute modes coefficients
[xMode,yMode] = fromGridToModes(x,y,PARAMnewton);

%compute grid points on the new grid
if PARAMnewton.legendre==1||PARAMnewton.legendre==2
    
    x = interpLegendreZeroOne(PARAM.t,xMode);
    y = interpLegendreZeroOne(PARAM.t,yMode);
    
elseif PARAMnewton.legendre==0
    
    x = chebcoeffs2chebvals(xMode);
    y = chebcoeffs2chebvals(yMode);
    
    x = chebfun(x,[0 1]);
    y = chebfun(y,[0 1]);
    
    x = x(PARAM.t);
    y = y(PARAM.t);
    
end












