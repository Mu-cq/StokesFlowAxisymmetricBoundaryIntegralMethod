%pload initial shape form continuation method

function [x,y,PARAMedge] = uploadFromEdgeTrackingBEM(PARAM)

%load data
fileLoad = [PARAM.res '/edgeTrackingDropBEM_edgeLoop=' num2str(PARAM.edgeLoopUpload) '_deltaEdge=' num2str(PARAM.deltaEdgeUpload) '_ODE=0_BC=' num2str(PARAM.BCupload) '_Bo=' num2str(PARAM.BondUpload) '_visc=' num2str(PARAM.viscUpload) '_n=' num2str(PARAM.elemUpload) '_maxDT=' num2str(PARAM.dtUpload) '.mat'];
edgeResults = load(fileLoad);
PARAMedge = edgeResults.PARAM;

%upload right shape
T = edgeResults.Tedge{PARAM.whichLoopUpload};
Y = edgeResults.Yedge{PARAM.whichLoopUpload};
VVV = edgeResults.Vedge{PARAM.whichLoopUpload};

%choose time with smallest residuals
res = zeros(numel(T),1);
for i = 1:numel(T)
    if PARAMedge.dropFrame==1
        Vhere = VVV{i};
    else
        error('Not implemented')
    end
    UxN = Vhere(1:2:end-1);
    UyN = Vhere(2:2:end);
    UnABS = sqrt(UxN.^2+UyN.^2);
    res(i) = norm(UnABS,Inf);
end

[minRes,ind] = min(res);

disp(['Shape is upload with res=' num2str(minRes) ' when base state is stable base state is below ' num2str(PARAMedge.resConverge)]);

Yhere = Y{ind};
x = Yhere(1:2:end-1)';
y = Yhere(2:2:end)';












