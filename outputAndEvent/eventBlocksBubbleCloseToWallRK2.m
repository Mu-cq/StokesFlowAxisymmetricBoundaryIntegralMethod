%activate spectral remesh

function value = eventBlocksBubbleCloseToWallRK2(t,var,T,Y,V,PARAM)

value = 0;
if isempty(var)==0
    
    %save results
    if PARAM.SaveDataIte==1
        
       %save results
       disp('Save results')
       cd(PARAM.res)
       save(PARAM.filename)
       cd(PARAM.here)
        
    end
    
    nNodes = (numel(var))/2;
    xBubble = var(1:2:2*nNodes-1);
    yBubble = var(2:2:2*nNodes);
    
    if min(yBubble)<0 || max(yBubble)>10 || min(xBubble)<PARAM.posWall
        disp('Interface has escaped the domain')
        value = 1;
    end
    
    if value==0

        fprintf('T=%f, gap=%1.6f \n',t,min(xBubble)-PARAM.posWall);

    end

end









