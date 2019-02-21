%activate spectral remesh

function value = eventBlocksODE(t,var,flag,tParametricBase,fCompute,beta,Hcc,cycleType,PARAM)

if isempty(var)==0
    
    %get position of the motor and bubble and bubble radius
    posMotor = var(1);
    posBubble = var(2);
    rBubble = var(3);
    
    %compute concentration along the axis
    [~,yLaplace,~,xVel,yVel,PARAMvel,Qbubble] = fCompute(t,var);
    Xline = linspace(-10+posMotor,20+posMotor,200);
    Yline = zeros(1,numel(Xline));
    [~,~,PHIfield] = computeConcentrationField(Xline,Yline,xVel,yVel,yLaplace,PARAMvel,0,1);
    
    %bubble concentration given by Henry's law
    phiBubble = Hcc*(beta+2/rBubble);

%     figure(1)
%     plot(Xline,PHIfield)
%     xlabel('z')
%     ylabel('c')
%     grid on
%     axis([min(Xline) max(Xline) 0 25])
%     drawnow
    
    %get right parameters
    PARAM.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
    PARAM.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
    PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;
    PARAM.x0_Circle(5) = posBubble;
    PARAM.rArc(5) = rBubble;

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

    if value==0

        fprintf('T=%f c_max=%f \n',t(1),max(PHIfield));

    end

    value = (value>0);

    if cycleType==1
        
        xcm = (max(x{5})+min(x{5}))/2;
        if xcm>max(x{4})
           disp('Bubble has exited the microrocket, simulation stops')
           value = 1;
        end
        
    elseif cycleType==2 || cycleType==3
    
        if max(PHIfield)>PARAM.maxConcAxis

            [maxPHI,indPHI] = max(PHIfield);

            z0 = Xline(indPHI);

            if abs(z0-posBubble-rBubble)>PARAM.distCritNucleation
                display(['New bubble nucleates in z=' num2str(z0) ' where c=' num2str(maxPHI)])
                value = 1;
            end

        end
    
    end

else
    value = 0;
end









