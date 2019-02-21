%activate spectral remesh

function value = eventBlocksTwoBubblesODE(t,var,flag,tParametricBase,fCompute2bubbles,endCycle,PARAM)

if isempty(var)==0
    
    %get position of the motor and bubble and bubble radius
    posMotor = var(1);
    posBubble1 = var(2);
    posBubble2 = var(3);
    rBubble1 = var(4);
    rBubble2 = var(5);
    
    %compute concentration along the axis
    [~,~,yLaplace,Qmass5,Qmass6,xVel,yVel,PARAMvel] = fCompute2bubbles(t,var);
    Xline = linspace(-10+posMotor,20+posMotor,200);
    Yline = zeros(1,numel(Xline));
    [~,~,PHIfield] = computeConcentrationField(Xline,Yline,xVel,yVel,yLaplace,PARAMvel,0,1);
    
    %get right parameters
    PARAM.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
    PARAM.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
    PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;
    PARAM.x0_Circle(5) = posBubble1;
    PARAM.rArc(5) = rBubble1;
    PARAM.x0_Circle(6) = posBubble2;
    PARAM.rArc(6) = rBubble2;

    %build geometry
    [x,y] = buildGeometryPanelsParametric(tParametricBase,PARAM);

    %check compenetration
    PARAM.typeBC = PARAM.typeBCstokes;
    PARAM.orderVariable = PARAM.orderVariableStokes;
    PARAM.orderGeometry = PARAM.orderGeometryStokes;
    value = 0;
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
    
    %stop if new bubble nucleates
    if max(PHIfield)>PARAM.maxConcAxis
        
        [maxPHI,indPHI] = max(PHIfield);
        
        z0 = Xline(indPHI);
        
        %closest bubble
        minDist = min([abs(z0-posBubble1-rBubble1) abs(z0-posBubble2-rBubble2)]);
        
        if minDist>PARAM.distCritNucleation
            display(['New bubble nucleates in z=' num2str(z0) ' where c=' num2str(maxPHI)])
            value = 1;
        end
        
    end
    
    %stop if larger bubble shrinks
    if Qmass5<0 && endCycle==2
        
        disp('Larger bubble shrinks')
        value = 1;
        
    end
    
    %stop if small bubble vanishes
    if Qmass6<0 && rBubble2<0.3
        
        disp('Smaller bubble shrinks and vanishes')
        value = 1;
        
    end
    
    if value==0

        fprintf('T=%f c_max=%f Q1=%f Q2=%f \n',t,max(PHIfield),Qmass5,Qmass6);

    end

else
    value = 0;
end









