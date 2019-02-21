%remesh panels with distribution

function [value,isterminal,direction] = activateRemeshDistrPanels(t,var,PARAM,panel,TTT)

direction = 0;
isterminal = 1;
value = 0;

%this is ook for one bubble with moving wall
nNodes = floor(numel(var)/2);
xDrop = var(1:2:2*nNodes-1)';
yDrop = var(2:2:2*nNodes)';

%kind of distribution
if PARAM.distr(panel)==0
    %uniform distribution
elseif PARAM.distr(panel)==2
    %distance from wall
else
    error('Not implemented')
end

%compute arc lenght
dx = diff(xDrop);   dy = diff(yDrop);
dl = sqrt(dx.^2+dy.^2);
dlMax = max(dl);    dlMin = min(dl);

if dlMax>PARAM.maxElem(panel) && min(abs(TTT-t))>(TTT(2)-TTT(1))*0.2 && size(var,2)==1
    
    value = 1;
    
elseif dlMin<PARAM.minSizeElemRemesh(panel) && min(abs(TTT-t))>(TTT(2)-TTT(1))*0.2 && size(var,2)==1
    
    value = 1;
    
end