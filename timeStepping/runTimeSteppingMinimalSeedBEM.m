%chose time stepping algorithm

function AreaAtTimeHorizon = runTimeSteppingMinimalSeedBEM(initialXY,Tsave,dt,initialDT,PARAM,V0)

%function for drop velocity
fOneDrop = @(t,var) computeVelocityDrop(t,var,PARAM);

%remesh function
remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);

%volume correction function
volCorrFunction = @(t,var) volCorrBubble(t,var,V0);

figure(99)
plot(initialXY(1:2:end-1),initialXY(2:2:end))
hold on
axis equal
grid on    
drawnow

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
    [~,Y] = RK2mostGeneralBubbleVolCorr(fOneDrop,Tsave,initialXY,dt,initialDT,outFun,remeshFunction,volCorrFunction);
elseif PARAM.ODE~=0
    [~,Y] = odeMatlabRemeshVolCorrBubble(fOneDrop,Tsave,initialXY,outFun,remeshFunction,volCorrFunction,options,ODE,PARAM);
end

%final shape and its area
Y = Y{end};
x = Y(1:2:end-1);   y = Y(2:2:end);
AreaAtTimeHorizon = -surf_gauss_vect(x',y');