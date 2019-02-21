%pload initial shape form continuation method

function [x,y,xBase,yBase,perturb,PARAMdns] = uploadFromDNS(PARAM)

%load data
fileLoad = [PARAM.res '/DropSpectral_ODE=' num2str(PARAM.ODE) '_Legendre=' num2str(PARAM.legendreUpload) '_BC=' num2str(PARAM.BCupload) '_Ca=' num2str(PARAM.CaUpload) '_visc=' num2str(PARAM.viscUpload) '_n=' num2str(PARAM.elemUpload) '_D=' num2str(PARAM.Dupload) '_maxDT=' num2str(PARAM.dtUpload) '_VolCorr=' num2str(PARAM.volCorrUpload) '_Tend=' num2str(PARAM.TendUpload) '.mat'];
dnsResults = load(fileLoad);
PARAMdns = dnsResults.PARAM;

%upload right shape
T = dnsResults.T;
Y = dnsResults.Y;
[minCheck,ind] = min(abs(T-PARAM.Tupload));
xMode = Y(ind,1:2:end-1)';
yMode = Y(ind,2:2:end)';
[x,y] = fromModesToGrid(xMode,yMode,PARAMdns);

%check if requested time is diiferent from uploaded time
if minCheck>1e-4
    error('No data at the requested time in the previous simulation')
end

%modes and perturbation of the base shape, this is useful for stability analysis
xBase = xMode;
yBase = yMode;
perturb = zeros(numel(PARAM.dealiasing)-1,1);

%stability works well only with the same spectral grid and same number of
%harmonics
if PARAM.algorithm==5 && (PARAM.legendre~=PARAMdns.legendre)
    error('Spectral basis must be the same as for the computation of the base state')
elseif PARAM.algorithm==5 && (PARAM.dealiasing~=PARAMdns.dealiasing)
    error('Number of harmonics must be the same as in the base shape')
end

%compute grid points on the new grid
if PARAMdns.legendre==1||PARAMdns.legendre==2
    
    x = interpLegendreZeroOne(PARAM.t,xMode);
    y = interpLegendreZeroOne(PARAM.t,yMode);
    
elseif PARAMdns.legendre==0
    
    x = chebcoeffs2chebvals(xMode);
    y = chebcoeffs2chebvals(yMode);
    
    x = chebfun(x,[0 1]);
    y = chebfun(y,[0 1]);
    
    x = x(PARAM.t);
    y = y(PARAM.t);
    
end












