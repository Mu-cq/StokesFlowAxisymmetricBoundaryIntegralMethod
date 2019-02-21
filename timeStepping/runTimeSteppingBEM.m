%chose time stepping algorithm

function allRes = runTimeSteppingBEM(Tsave,dt,initialDT,PARAM,V0)

%print to screen
printToScreenOneDrop(PARAM.Ca,PARAM.Bond,PARAM.visc(1),PARAM.dropFrame,PARAM.algorithm,[],[])

%initial condition
[initialXY,PARAM] = initialConditionDropBEM(PARAM);

%function for drop velocity
fOneDrop = @(t,var) computeVelocityDrop(t,var,PARAM);

%remesh function
remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);

%volume correction function
volCorrFunction = @(t,var) volCorrBubble(t,var,V0);

%output functions and events
if PARAM.ODE~=0
    outFun = @(t,var,flag) eventBlocksOneDrop(t,var,flag,fOneDrop,PARAM);
    %eventRemeshVolCorr = @(t,var) activateRemeshVolCorrDistrPanels(t,var,PARAM,4,Tsave,4/3*pi*alpha^3,volTol);
    eventRemeshVolCorr = @(t,var) activateRemeshVolCorrProximity(t,var,PARAM,4,Tsave,4/3*pi*alpha^3,volTol);
    options = odeset('RelTol',tol,'AbsTol',1e-3*tol,'OutputFcn',outFun,'Stats','off','MaxStep',dt,'Events',eventRemeshVolCorr,'InitialStep',initialDT);
else
    %output, saving and event function
    outFun = @(t,var,T,Y,V) eventBlocksOneDropRK2(t,var,T,Y,V,PARAM);
end

%time stepping
if PARAM.ODE==0
    [T,Y,V] = RK2mostGeneralBubbleVolCorr(fOneDrop,Tsave,initialXY,dt,initialDT,outFun,remeshFunction,volCorrFunction);
    allRes = {T Y V};
elseif PARAM.ODE~=0
    [T,Y] = odeMatlabRemeshVolCorrBubble(fOneDrop,Tsave,initialXY,outFun,remeshFunction,volCorrFunction,options,ODE,PARAM);
    allRes = {T Y};
end