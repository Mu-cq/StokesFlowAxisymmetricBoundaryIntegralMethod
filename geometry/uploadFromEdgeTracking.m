%pload initial shape form continuation method

function [x,y,xBase,yBase,perturb,PARAMedge] = uploadFromEdgeTracking(PARAM)

%load data
fileLoad = [PARAM.res '/edgeTrackingDrop_edgeLoop=' num2str(PARAM.edgeLoopUpload) '_deltaEdge=' num2str(PARAM.deltaEdgeUpload) '_ODE=2_Legendre=' num2str(PARAM.legendreUpload) '_BC=' num2str(PARAM.BCupload) '_Ca=' num2str(PARAM.CaUpload) '_visc=1_n=' num2str(PARAM.elemUpload) '_maxDT=' num2str(PARAM.dtUpload) '_VolCorr=' num2str(PARAM.volCorrUpload) '.mat'];
edgeResults = load(fileLoad);
PARAMedge = edgeResults.PARAM;

%upload right shape
T = edgeResults.Tedge{PARAM.whichLoopUpload};
Y = edgeResults.Yedge{PARAM.whichLoopUpload};
[~,ind] = min(abs(T-PARAM.Tupload));
xMode = Y(ind,1:2:end-1)';
yMode = Y(ind,2:2:end)';
[x,y] = fromModesToGrid(xMode,yMode,PARAMedge);

%modes and perturbation of the base shape, this is useful for stability analysis
xBase = xMode;
yBase = yMode;
perturb = zeros(numel(PARAM.dealiasing)-1,1);

%stability works well only with the same spectral grid and same number of
%harmonics
if PARAM.algorithm==5 && (PARAM.legendre~=PARAMedge.legendre)
    error('Spectral basis must be the same as for the computation of the base state')
elseif PARAM.algorithm==5 && (PARAM.dealiasing~=PARAMedge.dealiasing)
    error('Number of harmonics must be the same as in the base shape')
end

%compute grid points on the new grid
if PARAMedge.legendre==1||PARAMedge.legendre==2
    
    x = interpLegendreZeroOne(PARAM.t,xMode);
    y = interpLegendreZeroOne(PARAM.t,yMode);
    
elseif PARAMedge.legendre==0
    
    x = chebcoeffs2chebvals(xMode);
    y = chebcoeffs2chebvals(yMode);
    
    x = chebfun(x,[0 1]);
    y = chebfun(y,[0 1]);
    
    x = x(PARAM.t);
    y = y(PARAM.t);
    
end












