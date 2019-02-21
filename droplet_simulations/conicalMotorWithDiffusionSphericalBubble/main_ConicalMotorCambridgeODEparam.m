%compute the motion of a conical motor with diffusion of chemical species

function main_ConicalMotorCambridgeODEparam(dt,tol,ODE,thetaCone,L,repOn,repIntensity,BN,Hcc,beta,cycleType,Tend,checkpoint,results)

%results
PARAM.res = results;
PARAM.here = pwd;

%physical parameters
r = 1;                      % average cone radius
h = 0.2;                    % wall thickness

%nucleation
bubbleNucleates = BN;
PARAM.distCritNucleation = 2;
PARAM.smallestRadius = 0;
PARAM.maxRadiusNucleation = 1;

%solve Laplace
PARAM.SolveLaplace = 1;

%options
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.STlaplace = 1;

%time discretization
Tstart = 0;

%saving options
SaveHowMany = checkpoint;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time

%geometry parameters
nPerLenght = 5;
PARAM.panels = [4 1];                              % panels per block
PARAM.rotate = [thetaCone thetaCone thetaCone thetaCone 0];
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) nan];                  % number of element per panel
PARAM.geometryPanel = [0 1 0 1 1 1];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi pi];               % theta starting point for the arc
PARAM.rArc = [nan h/2 nan h/2 nan nan];               % theta starting point for the arc
PARAM.x0_Circle = [nan 0 nan L nan nan];
PARAM.y0_Circle = [nan r nan r 0 0];
PARAM.xCrotate = [0 0 0 0 0 0];
PARAM.yCrotate = [r r r r 0 0];

%numerics parameters for Laplace
PARAM.typeBClaplace = [2 2 2 2 1 1];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0 0 0];     % 0 is straight, 1 is curved (spline)

%numerics parameters for Stokes
PARAM.typeBCstokes = [1 1 1 1 5];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 0 0];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [1 1 1 1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0 1];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [1 1];                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [1 1];
PARAM.deflationConstant = [4*pi 4*pi];
PARAM.addFlow = 0;

%repulsive forces
PARAM.repulsiveForces = [1 1];
PARAM.repulsiveOn = repOn;
PARAM.debeyLenght = 45;
PARAM.coeffRepulsive = repIntensity*sin(thetaCone);

%remesh
PARAM.remeshType = [0 0 1 1 1 1];         % 1 element is split into two parts when coming close to another block, 2 uses a distribution
PARAM.remeshStep = 2;
PARAM.remeshProximity = {[] [] 5 5 3 3};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = [0 0 L/PARAM.n(3) pi*h/PARAM.n(4) nan nan];
PARAM.distActivateRemesh = [1 1 1 1 1 1];
PARAM.adaptCoeff = [4 4 4 4 4 4];
PARAM.minSizeElemRemesh = [1e-3 1e-3 1e-3 1e-3 1e-3 1e-3]/2;
PARAM.coeffDist = 2;    % how much smaller than distnce it has to be

%function profile for BCs
PARAM.fluxBC = {0 0 -1 0 [] []};
PARAM.velBC = {0 0 0 0 []};
PARAM.stressBC{5} = 0;

%choose bubble position and size
[PARAM,valueNuc] = chooseBubblePositionAndSize(PARAM,beta,Hcc,nPerLenght,bubbleNucleates);

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)
printToScreenTimeStepping(ODE,dt,tol)
printToScreenRemesh(PARAM)
printToScreenConicalMotorWithSphere(thetaCone,r,L,PARAM.x0_Circle(5),PARAM.rArc(5),Hcc,beta)

%build geometry
[~,~,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

%initial condition
initial = [0 PARAM.x0_Circle(5) PARAM.rArc(5)]';

%create filename
PARAM.filename = chooseFilenameConicalMotor(PARAM,nPerLenght,beta,Hcc,bubbleNucleates,L,Tend,thetaCone,dt,ODE,tol,0,cycleType);

%function for motor velocity
fMotor = @(t,var) computeVelocityConicalMotorODE(t,var,tParametricBase,nPerLenght,PARAM,Hcc,beta);

%output and event function
outFun = @(t,var,flag) eventBlocksODE(t,var,flag,tParametricBase,fMotor,beta,Hcc,cycleType,PARAM);

%ODE options
options = odeset('RelTol',tol,'AbsTol',1e-6,'OutputFcn',outFun,'Stats','off','MaxStep',dt);

tic

if valueNuc==0

    %time stepping
    if ODE==0
        [T,Y,V] = RK2mostGeneral(fMotor,Tsave,initial,dt,dt,outFun);
    elseif ODE==1
        [T,Y] = ode45(fMotor,Tsave, initial, options);
    elseif ODE==2
        [T,Y] = ode23t(fMotor,Tsave, initial, options);
    elseif ODE==3
        [T,Y] = ode23s(fMotor,Tsave, initial, options);
    end

elseif valueNuc==1
    
    display('Bubble does not nucleate, it is too large or it is already shrinking')
    
end

simulationTime = toc;

clear fMotor
clear outFun

%save results
display('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

display('The end')


