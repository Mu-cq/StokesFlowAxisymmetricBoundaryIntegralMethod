%chose time stepping algorithm

function [Tedge,Yedge,Vedge] = edgeTrackingBEM(Tsave,dt,initialDT,PARAM,V0)

%loops
loops = PARAM.edgeLoop;                     % loop of edge tracking
loopStart = PARAM.edgeStartingLoop;         % starting loop of edge tracking

%initialize
Tedge = cell(loops,1);
Yedge = cell(loops,1);
Vedge = cell(loops,1);

%load data form previous simulation
if loopStart>1
    
    here = pwd;
    cd(PARAM.res)
    Told = load(PARAM.filename,'Tedge');
    Yold = load(PARAM.filename,'Yedge');
    Vold = load(PARAM.filename,'Vedge');
    cd(here)
    
    Tedge(1:loopStart-1) = Told.Tedge(1:loopStart-1);
    Yedge(1:loopStart-1) = Yold.Yedge(1:loopStart-1);
    Vedge(1:loopStart-1) = Vold.Vedge(1:loopStart-1);
    
end

for l = loopStart:loops

    %print to screen
    PARAM.loopNum = l;
    printToScreenOneDrop(PARAM.Ca,PARAM.Bond,PARAM.visc,PARAM.dropFrame,PARAM.algorithm,l,loops);

    %initial condition
    if l==1 || l==2
        %ellipse
        PARAM.D = PARAM.Dedge(l);
        [~,~,Lside] = ellipseCartesian(linspace(0,pi,100),PARAM.D);
        PARAM.maxElem = 1.1*Lside/PARAM.n(1);
        initialXY = initialConditionDropBEM(PARAM);
        TsaveNow = Tsave;
    elseif l>2
        
        if PARAM.interpOrExtrap==1
            
            %build shape from previous loops by interpolation
            [Tinitial,initialXY] = buildShapeEdgeBisectionBEM(Tedge(1:l-1),Yedge(1:l-1),PARAM,V0);
        
        elseif PARAM.interpOrExtrap==2
            
            %build shape from previous loops by extrapolation
            [Tinitial,initialXY] = buildShapeEdgeExtrapolateBEM(Tedge(1:l-1),Yedge(1:l-1),PARAM,V0);
            
        else

            error('Not implemented')
        
        end
        
        %time for next simulation
        TsaveNow = Tinitial+Tsave;
        
    end

    %function for drop velocity
    fOneDrop = @(t,var) computeVelocityDrop(t,var,PARAM);

    %remesh function
    remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);

    %volume correction function
    volCorrFunction = @(t,var) volCorrBubble(t,var,V0);

    %output functions and events
    if PARAM.ODE~=0
        outFun = @(t,var,flag) eventBlocksOneDrop(t,var,flag,fOneDrop,PARAM);
        eventRemeshVolCorr = @(t,var) activateRemeshVolCorrProximity(t,var,PARAM,4,Tsave,4/3*pi*alpha^3,volTol);
        options = odeset('RelTol',PARAM.tol,'AbsTol',1e-3*PARAM.tol,'OutputFcn',outFun,'Stats','off','MaxStep',dt,'Events',eventRemeshVolCorr,'InitialStep',initialDT);
    else
        %output, saving and event function
        outFun = @(t,var,T,Y,V) eventBlocksOneDropRK2(t,var,T,Y,V,PARAM);
    end

    %time stepping
    if PARAM.ODE==0
        [T,Y,V] = RK2mostGeneralBubbleVolCorr(fOneDrop,TsaveNow,initialXY,dt,initialDT,outFun,remeshFunction,volCorrFunction);
    elseif PARAM.ODE~=0
        [T,Y] = odeMatlabRemeshVolCorrBubble(fOneDrop,TsaveNow,initialXY,outFun,remeshFunction,volCorrFunction,options,PARAM.ODE,PARAM);
    end
    
    Tedge{l} = T;
    Yedge{l} = Y;
    Vedge{l} = V;
    
    %save results
    if PARAM.saveData==1
        here = pwd;
        cd(PARAM.res)
        disp('Save data')
        save(PARAM.filename)
        cd(here)
    end

end

