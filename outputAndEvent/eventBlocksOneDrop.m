%activate spectral remesh

function value = eventBlocksOneDrop(t,var,flag,fCompute,PARAM)

value = 0;
if isempty(var)==0 && numel(t)==1
    
    nNodes = ceil((numel(var)-1)/2);
    xBubble = var(1:2:2*nNodes-1);
    yBubble = var(2:2:2*nNodes);
    PARAM.n(1) = numel(xBubble)-1;
    
    if min(yBubble)<-1e-10 || max(yBubble)>10
        disp('Interface has escaped the domain')
        value = 1;
    end
    
    if isempty(PARAM.Ca) && max(xBubble)-min(xBubble)>3
        disp('tail is elongating in an unstable manner')
        value = 1;
    end
    
    %check residuals
    res = 1;
    if t(1)>0 && value==0
        Vhere = fCompute(t,var);
        res = NormalVelocityDropFrame(xBubble',yBubble',Vhere(1:2:end-1),Vhere(2:2:end));
        res = norm(res,Inf);
    end
    
    if value==0

        fprintf('T=%f, res=%1.10f \n',t,res);

    end

end









