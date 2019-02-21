%pload initial shape form continuation method

function [x,y,xBase,yBase,perturb,PARAMnewton] = uploadFromNewton(PARAM)

%load data
here = pwd;
cd(PARAM.res)
fileLoad = ['newtonMethodSpectralXYmodes_n=' num2str(PARAM.elemUpload) '_Ca=' num2str(PARAM.CaUpload) '_visc=' num2str(PARAM.viscUpload) '_Legendre=' num2str(PARAM.legendreUpload) '_BC=' num2str(PARAM.BCupload) '.mat'];
newtonResults = load(fileLoad);
cd(here)
PARAMnewton = newtonResults.PARAM;

%check mappings
printError = 0;
try
    if (PARAM.chooseRemeshMapping~=PARAMnewton.chooseRemeshMapping) || (PARAM.howStrong~=PARAMnewton.howStrong)
        printError = 1;
    end
catch
    disp('No mapping was applied on the updated shape or this "check" was not implemented yet. Be careful not to use a different mapping in the new simulation')
end
if printError==1
    error('Mesh mapping has to be the same otherwise the shape will be deformed')
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

%compute modes from grid
x = newtonResults.x;
y = newtonResults.y;
[xMode,yMode] = fromGridToModes(x,y,PARAMnewton);

%compute grid points on the new grid
if PARAMnewton.legendre==1||PARAMnewton.legendre==2
    
    x = interpLegendreZeroOne(PARAM.t,xMode);
    y = interpLegendreZeroOne(PARAM.t,yMode);
    
elseif PARAM.legendre==0
    
    x = chebcoeffs2chebvals(xMode);
    y = chebcoeffs2chebvals(yMode);
    
    x = chebfun(x,[0 1]);
    y = chebfun(y,[0 1]);
    
    x = x(PARAM.t);
    y = y(PARAM.t);
    
end












