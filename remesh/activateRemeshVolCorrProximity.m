%remesh panels with distribution

function [value,isterminal,direction] = activateRemeshVolCorrProximity(t,var,PARAM,panel,TTT,Vin,tolVol)

direction = 0;
isterminal = 1;
value = 0;

%this is ook for one bubble with moving wall
%         
nNodes = floor(numel(var)/2);
xDrop = var(1:2:2*nNodes-1)';
yDrop = var(2:2:2*nNodes)';

%check if to activate remesh
[~,~,checkActiveRemesh] = remeshBubbleSplitPanels([],xDrop,yDrop,PARAM,panel);

if checkActiveRemesh>=1 && min(abs(TTT-t))>(TTT(2)-TTT(1))*0.2 && size(var,2)==1
    
    value = 1;
    
end

%check volume
Vnow = axis_int_gauss_vect(xDrop,yDrop);
errV = abs(Vnow-Vin)/Vin;
if errV>tolVol
    
    %display('Volume correction')
    value = 1;
    
end