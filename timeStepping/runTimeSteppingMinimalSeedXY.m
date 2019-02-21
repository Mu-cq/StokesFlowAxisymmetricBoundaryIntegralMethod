%chose time stepping algorithm

function [Aperturb,T,Y] = runTimeSteppingMinimalSeedXY(xy,Tsave,maxDT,initialDT,PARAM,V0,Abase)

%initial condition
xGrid = xy(1:2:end-1);  yGrid = xy(2:2:end);
[xMode,yMode] = fromGridToModes(xGrid,yGrid,PARAM);
initial = zeros(2*numel(xMode),1);
initial(1:2:end-1) = xMode;
initial(2:2:end) = yMode;

figure(100)
plot(xGrid,yGrid)
hold on
axis equal
drawnow

%output and event function
myOut1 = @(t,y,flag) checkCoeffsAndConvergence(t,y,flag,PARAM,V0);
eventRemesh = @(t,y) eventSpectralRemesh(t,y,PARAM);

%nonlinear function
f = functionDropSpectral(PARAM);

%ODE options
options = odeset('RelTol',PARAM.tresholdCoeffs*PARAM.Rescale,'AbsTol',PARAM.tresholdCoeffs*PARAM.Rescale,'OutputFcn',myOut1,'Stats','off','MaxStep',maxDT,'Events',eventRemesh);

if PARAM.ODE==1
    if PARAM.remesh==0
        [T,Y] = ode45(f,Tsave, initial, options);
    elseif PARAM.remesh==1
        [T,Y] = remeshODE45matlab(f,Tsave,initial,options,PARAM);
    end
elseif PARAM.ODE==2
    [T,Y] = RK2(f,Tsave,initial,maxDT,initialDT,myOut1,PARAM);
elseif PARAM.ODE==3
    [T,Y] = ode23s(f,Tsave, initial, options);
elseif PARAM.ODE==4
    [T,Y] = ode23(f,Tsave, initial, options);
elseif PARAM.ODE==5
    [T,Y] = ode113(f,Tsave, initial, options);
elseif PARAM.ODE==6
    if PARAM.remesh==0
        [T,Y] = ode23t(f,Tsave, initial, options);
    elseif PARAM.remesh==1
        [T,Y] = remeshODE23tmatlab(f,Tsave,initial,options,PARAM);
    end
elseif PARAM.ODE==7
    [T,Y] = ode15s(f,Tsave, initial, options);
elseif PARAM.ODE==8
    [T,Y] = ode23tb(f,Tsave, initial, options);
elseif PARAM.ODE==9
    [T,Y] = RK4(f,Tsave,initial,maxDT,initialDT,myOut1,PARAM);
end

%compute surface area perturbation (compared to base state) at last iteration
xMode = Y(end,1:2:end-1)';
yMode = Y(end,2:2:end)';
[xGrid,yGrid] = fromModesToGrid(xMode,yMode,PARAM);
Aperturb = -surfaceCurvilinearAxisSpectral(xGrid,yGrid,PARAM)+Abase;




