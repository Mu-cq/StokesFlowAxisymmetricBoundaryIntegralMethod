%activate spectral remesh

function value = eventBlocksOneDropRK2(t,var,T,Y,V,PARAM)

value = 0;
if isempty(var)==0
    
    xBubble = var(1:2:end-1);
    yBubble = var(2:2:end);
    
    %save results
    if PARAM.SaveDataIte==1
        
       %save results
       disp('Save results')
       cd(PARAM.res)
       save(PARAM.filename)
       cd(PARAM.here)
        
    end
    
    if min(yBubble)<0 || max(yBubble)>10
        disp('Interface has escaped the domain')
        value = 1;
    end
    
    %check residuals
    res = 1;
    if t(1)>0 && value==0
        tHere = t(1);
        [~,ind] = min(abs(T-tHere));
        Vhere = V{ind};
        Yhere = Y{ind}';
        res = NormalVelocityDropFrame(Yhere(1:2:end-1),Yhere(2:2:end),Vhere(1:2:end-1),Vhere(2:2:end));
        res = norm(res,Inf);
    end
    
    if isempty(PARAM.Ca) && res>1.5 && t>0
        disp('Instabilities are developing')
        value = 1;
    end
    
    if res<PARAM.resConverge
        disp('Convergence has been reached')
        value = 1;
    end
    
    if value==0

        fprintf('T=%f, res=%1.10f \n',t,res);

    end

end









