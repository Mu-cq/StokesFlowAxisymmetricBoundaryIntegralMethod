%activate spectral remesh

function value = eventBlocksBubbleCloseToWall(t,var,flag,PARAM)

value = 0;
if isempty(var)==0 && numel(t)==1
    
    nNodes = ceil((numel(var)-1)/2);
    xBubble = var(1:2:2*nNodes-1);
    yBubble = var(2:2:2*nNodes);
    PARAM.n(4) = numel(xBubble)-1;
    
    if min(yBubble)<0 || max(yBubble)>10 || min(xBubble)<PARAM.posWall
        disp('Interface has escaped the domain')
        value = 1;
    end
    
    if value==0

        fprintf('T=%f, gap=%1.6f \n',t,min(xBubble)-PARAM.posWall);

    end

end









