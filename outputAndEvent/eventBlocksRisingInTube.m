%activate spectral remesh

function value = eventBlocksRisingInTube(t,var,flag,tParametricBase,fCompute,PARAM)

value = 0;
if isempty(var)==0 && numel(t)==1
    
    %save results
%     if PARAM.SaveDataIte==1
%         
%        %save results
%        disp('Save results')
%        cd(PARAM.res)
%        save(PARAM.filename)
%        cd(PARAM.here)
%         
%     end
    
    nNodes = ceil((numel(var)-1)/2);
    xBubble = var(1:2:2*nNodes-1);
    yBubble = var(2:2:2*nNodes);
    PARAM.n(4) = numel(xBubble)-1;
    
    if min(yBubble)<-1e-10 || max(yBubble)>10
        disp('Interface has escaped the domain')
        value = 1;
    end

    %build geometry
    [x,y] = buildGeometryPanelsParametric(tParametricBase,PARAM);
    x{4} = xBubble';
    y{4} = yBubble';

    %check compenetration
    PARAM.typeBC = PARAM.typeBCstokes;
    PARAM.orderVariable = PARAM.orderVariableStokes;
    PARAM.orderGeometry = PARAM.orderGeometryStokes;
    for i = 1:numel(PARAM.blockType)

        [xHere1,yHere1] = getBlockCoordinates(x,y,PARAM,i);

        for k = 1:numel(PARAM.blockType)

            if k~=i
                [xHere2,yHere2] = getBlockCoordinates(x,y,PARAM,k);
                InOut = FigInOut(xHere1,yHere1,xHere2',yHere2');
                if max(InOut)>pi
                    warning('Compenetration occurs')
                    value = value+1;
                end
            end

        end

    end
    value = (value>0);
    
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









