%chose time stepping algorithm

function [Tedge,Yedge] = edgeTracking(Tsave,maxDT,initialDT,PARAM,V0)

%loops
loops = PARAM.edgeLoop;                     % loop of edge tracking
loopStart = PARAM.edgeStartingLoop;         % starting loop of edge tracking

%initialize
Tedge = cell(loops,1);
Yedge = cell(loops,1);

%load data form previous simulation
if loopStart>1
    
    here = pwd;
    cd(PARAM.res)
    Told = load(PARAM.filename,'Tedge');
    Yold = load(PARAM.filename,'Yedge');
    cd(here)
    
    Tedge(1:loopStart-1) = Told.Tedge(1:loopStart-1);
    Yedge(1:loopStart-1) = Yold.Yedge(1:loopStart-1);
    
end

for l = loopStart:loops

    %print to screen
    PARAM.loopNum = l;
    printToScreen(PARAM);

    %initial condition
    if l==1 || l==2
        %ellipse
        PARAM.D = PARAM.Dedge(l);
        PARAM.overlapMode = PARAM.overlapModeEdge(l);
        PARAM.whichmode = PARAM.whichModeEdge(l);
        initial = initialConditionDrop(PARAM);
        TsaveNow = Tsave;
    elseif l>2
        %build shape from previous loops
        [Tinitial,initial] = buildShapeEdgeBisection(Tedge(1:l-1),Yedge(1:l-1),PARAM,V0);
        
        %time for next simulation
        TsaveNow = Tinitial+Tsave;
        
    end

    %output and event function
    myOut1 = @(t,y,flag) checkCoeffsAndConvergence(t,y,flag,PARAM,V0);
    eventRemesh = @(t,y) eventSpectralRemesh(t,y,PARAM);

    %nonlinear function
    f = functionDropSpectral(PARAM);

    %ODE options
    options = odeset('RelTol',PARAM.tresholdCoeffs*PARAM.Rescale,'AbsTol',PARAM.tresholdCoeffs*PARAM.Rescale,'OutputFcn',myOut1,'Stats','off','MaxStep',maxDT,'Events',eventRemesh);

    %solve ode fro time marching
    if PARAM.ODE==1
        if PARAM.remesh==0
            [T,Y] = ode45(f,TsaveNow, initial, options);
        elseif PARAM.remesh==1
            [T,Y] = remeshODE45matlab(f,TsaveNow,initial,options,PARAM);
        end
    elseif PARAM.ODE==2
        [T,Y] = RK2(f,TsaveNow,initial,maxDT,initialDT,myOut1,PARAM);
    elseif PARAM.ODE==3
        [T,Y] = ode23s(f,TsaveNow, initial, options);
    elseif PARAM.ODE==4
        [T,Y] = ode23(f,TsaveNow, initial, options);
    elseif PARAM.ODE==5
        [T,Y] = ode113(f,TsaveNow, initial, options);
    elseif PARAM.ODE==6
        if PARAM.remesh==0
            [T,Y] = ode23t(f,TsaveNow, initial, options);
        elseif PARAM.remesh==1
            [T,Y] = remeshODE23tmatlab(f,TsaveNow,initial,options,PARAM);
        end
    elseif PARAM.ODE==7
        [T,Y] = ode15s(f,TsaveNow, initial, options);
    elseif PARAM.ODE==8
        [T,Y] = ode23tb(f,TsaveNow, initial, options);
    elseif PARAM.ODE==9
        [T,Y] = RK4(f,TsaveNow,initial,maxDT,initialDT,myOut1,PARAM);
    end
    
    Tedge{l} = T;
    Yedge{l} = Y;
    
    %save results
    if PARAM.saveData==1
        here = pwd;
        cd(PARAM.res)
        display('Save data')
        save(PARAM.filename)
        cd(here)
    end

end