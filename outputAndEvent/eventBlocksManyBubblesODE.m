%activate spectral remesh

function [value,newBubble,z0,r0,indOut] = eventBlocksManyBubblesODE(t,var,flag,tParametricBase,fCompute,beta,Hcc,PARAM)

newBubble = [];
r0 = [];
z0 = [];
indOut = [];
if isempty(var)==0
    
    nBubble = (numel(var)-1)/2;
    nWall = numel(PARAM.n)-nBubble;
    
    %get position of the motor and bubble and bubble radius
    posMotor = var(1);

    %for many bubbles
    posBubble = zeros(nBubble,1);
    rBubble = zeros(nBubble,1);
    for i = 1:nBubble
        posBubble(i) = var(i+1);
        rBubble(i) = var(1+nBubble+i);
        PARAM.x0_Circle(i+nWall) = posBubble(i);
        PARAM.rArc(i+nWall) = rBubble(i);
    end
    
    %compute concentration along the axis
    [~,~,yLaplace,Qmass,xVel,yVel,PARAMvel] = fCompute(t,var);
    Xline = linspace(-10+posMotor,20+posMotor,200);
    Yline = zeros(1,numel(Xline));
    [~,~,PHIfield] = computeConcentrationField(Xline,Yline,xVel,yVel,yLaplace,PARAMvel,0,1);
    
    %get right parameters
    PARAM.xStart(1:nWall) = PARAM.xStart(1:nWall)+posMotor;
    PARAM.xEnd(1:nWall) = PARAM.xEnd(1:nWall)+posMotor;
    PARAM.x0_Circle(1:nWall) = PARAM.x0_Circle(1:nWall)+posMotor;

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

        fprintf('T=%f c_max=%f \n',t,max(PHIfield));

    end

    value = (value>0);
    
    if max(PHIfield)>PARAM.maxConcAxis
        
        [maxPHI,indPHI] = max(PHIfield);
        
        z0 = Xline(indPHI);
        
        %closest bubble
        minDist = min(abs(z0-posBubble-rBubble));
        
        if minDist>PARAM.distCritNucleation
            display(['New bubble nucleates in z=' num2str(z0) ' where c=' num2str(maxPHI)])
            r0 = 2*Hcc/(maxPHI-beta*Hcc);
            r0 = r0+0.1;
            
            if r0<PARAM.smallestRadius
                
               r0 = PARAM.smallestRadius;
               display(['Minimum bubble size allowed is r=' num2str(PARAM.smallestRadius)])
                
            end
            
            newBubble = 1;
            value = 1;
        end
        
    end
    
    %stop if small bubble vanishes
    for i = 1:nBubble
        if Qmass(i)<0 && rBubble(i)<0.3

            disp('Small bubble shrinks and vanishes')
            newBubble = -1;
            indOut = i;
            value = 1;
            break

        end
    end

else
    value = 0;
end









