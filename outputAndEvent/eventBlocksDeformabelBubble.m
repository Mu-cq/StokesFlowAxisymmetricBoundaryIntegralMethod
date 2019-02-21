%activate spectral remesh

function value = eventBlocksDeformabelBubble(t,var,T,Y,V,tParametricBase,PARAM)

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
    
    nNodes = (numel(var)-1)/2;
    xBubble = var(1:2:2*nNodes-1);
    yBubble = var(2:2:2*nNodes);
    
    if min(yBubble)<0 || max(yBubble)>10
        disp('Interface has escaped the domain')
        value = 1;
    end
    
    %get position of the motor and bubble and bubble radius
    posMotor = var(end);
    
    %get right parameters
    PARAM.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
    PARAM.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
    PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;

    %build geometry
    [x,y] = buildGeometryPanelsParametric(tParametricBase,PARAM);
    x{5} = xBubble';
    y{5} = yBubble';

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
    
    %check is bubble is out and spherical
    aBubble = max(xBubble)-min(xBubble);
    bBubble = 2*max(yBubble);
    Dbubble = (aBubble-bBubble)/(aBubble+bBubble);
    if (max(xBubble)>max(x{4}) && Dbubble<2e-3) || (min(xBubble)<min(x{2}) && Dbubble<2e-3)
        disp('Bubble is out and spherical')
        value = 1;
    end
    
    if value>=1
        
        value = 1;
        
    end
    
    if value==0

        fprintf('T=%f, D=%f \n',t,Dbubble);

    end

end









