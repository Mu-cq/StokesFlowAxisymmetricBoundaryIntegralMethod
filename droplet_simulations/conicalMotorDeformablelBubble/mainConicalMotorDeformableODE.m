%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%results
PARAM.res = '../results';
PARAM.here = pwd;

%geometrical and physical parameters
r = 1;                      % radial coordinate for pivoting angle
L = 10;                     % motor lenght
h = 0.2;                    % wall thickness
thetaCone = 0;          % inclinitation of the cone
PARAM.visc(5) = 0;
PARAM.D(5) = 0;
Ca = 0.5;
alpha = 0.5;  squeezedBubble = 0;
beta = 1;
PARAM.massFlux = 1;
PARAM.ellipseShape = 1;
initialPosBubble = 1;    nDrop = 20;

%options
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];

%time discretization
Tstart = 0;
Tend = 2000;
dt = 5e-2;
ODE = 0;
tol = [];

%saving options
SaveHowMany = 1e3;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time
PARAM.SaveDataIte = 1;

%geometry parameters
nPerLenght = 5;
PARAM.panels = [4 1];                              % panels per block
PARAM.rotate = [thetaCone thetaCone thetaCone thetaCone 0]/180*pi;
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) nDrop];                  % number of element per panel
PARAM.geometryPanel = [0 1 0 1 2];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi];               % theta starting point for the arc
PARAM.rArc = [nan h/2 nan h/2 alpha];               % theta starting point for the arc
PARAM.x0_Circle = [nan 0 nan L initialPosBubble];
PARAM.y0_Circle = [nan r nan r 0];
PARAM.xCrotate = [0 0 0 0 0];
PARAM.yCrotate = [r r r r 0];

%numerics parameters for Stokes
PARAM.typeBCstokes = [0 0 0 0 2];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0 1];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0 1];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 0 2];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [1 1 1 1 2];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0 1];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [1 2];                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [1 1];
PARAM.deflationConstant = [4*pi 4*pi];
PARAM.addFlow = 0;
PARAM.blockVolCorr = [0 0];
%repulsive forces
PARAM.repulsiveForces = [4 4];
PARAM.repulsiveOn = 2e-1;
PARAM.coeffRepulsive = -10;
PARAM.smoothingRep = 1;

%remesh
PARAM.remeshType = [0 1 1 1 4];         % 1 element is split into two parts when coming close to another block, 2 uses a distribution, 4 remesh with a certain distributio (usually for deformable objects)
PARAM.remeshStep = 2;
PARAM.remeshProximity = {[] 5 5 5 3};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = [0 pi*h/PARAM.n(2) L/PARAM.n(3) pi*h/PARAM.n(4) alpha*pi/PARAM.n(5)];
PARAM.distActivateRemesh = [1 1 1 1 1 1];
PARAM.adaptCoeff = [4 4 4 4 4 4];
PARAM.minSizeElemRemesh = [1e-2 1e-2 1e-2 1e-2 1e-2];
PARAM.coeffDist = 2;    % how much smaller than distnce it has to be
PARAM.distr = [nan nan nan nan 0];
PARAM.maxNumberTotalElem = 1e3;

%function profile for BCs
PARAM.velBC = {0 0 0 0 0};
PARAM.velBCaxial = {0 0 0 0};
PARAM.velBCradial = {0 0 0 0};
PARAM.stressBC{5} = 1;

%print to screen
printToScreenStokes(PARAM)
printToScreenTimeStepping(ODE,dt,tol)
printToScreenRemesh(PARAM)
printToScreenConicalMotor(thetaCone,r,L,PARAM.x0_Circle(5),PARAM.rArc(5),Ca,PARAM.massFlux)

%build geometry
[xInitial,yInitial,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

%build squeezed droplet
if squeezedBubble==1
    
    [xInitial{5},yInitial{5}] = squeezedDropletFitLegendre(r,alpha,h,-2*h,thetaCone/180*pi,L/2,L,nDrop);
    xInitial{5} = xInitial{5}+L/2;
    
end

%initial condition
initialXY = zeros(2*numel(xInitial{5}),1);
initialXY(1:2:end-1) = xInitial{5};
initialXY(2:2:end) = yInitial{5};
initial = [initialXY; 0];

%create filename
PARAM.filename = chooseFilenameConicalMotorDeformable(PARAM,sum(PARAM.n),L,Tend,thetaCone,dt,PARAM.rArc(5),PARAM.x0_Circle(5),Ca,ODE,PARAM.massFlux);

%function for motor velocity
fMotor = @(t,var) computeVelocityDeformableBubbleODE(t,var,tParametricBase,Ca,beta,PARAM);

%remesh function
remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);

%output functions and events
if ODE~=0
    outFun = @(t,var,flag) eventBlocksDeformabelBubbleODE(t,var,flag,tParametricBase,PARAM);
    eventRemesh = @(t,var) activateRemeshDistrPanels(t,var,PARAM,5,Tsave);
    options = odeset('OutputFcn',outFun,'Stats','off','MaxStep',dt,'Events',eventRemesh,'InitialStep',dt*1e-2);
else
    %output, saving and event function
    outFun = @(t,var,T,Y,V) eventBlocksDeformabelBubble(t,var,T,Y,V,tParametricBase,PARAM);
end

tic

%time stepping
if ODE==0
    [T,Y,V] = RK2mostGeneralBubble(fMotor,Tsave,initial,dt,dt,outFun,remeshFunction);
elseif ODE~=0
    [T,Y] = odeMatlabRemeshBubble(fMotor,Tsave,initial,outFun,remeshFunction,options,ODE,PARAM);
end

simulationTime = toc;

%clear function handle
clear outfun
clear fMotor
clear remeshFunction
if ODE~=0
    clear eventRemesh
end

%save results
disp('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

disp('The end')


