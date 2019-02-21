%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/droplet_simulations/results';
PARAM.here = pwd;

%physical parameters
r = 1;                      % average cone radius
rBubble = 0.8;              % radius of the bubble
L = 10;                     % motor lenght
h = 0.2;                    % wall thickness
thetaCone = pi/64;          % inclinitation of the cone
xcmBubble = 5;
K = 0.375;
beta = 1.4;

%options
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.STlaplace = 1;
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];

%time discretization
Tstart = 0;
Tend = 0.4;
dt = 0.001;

%saving options
SaveHowMany = 200;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time

%geometry parameters
nPerLenght = 5;
PARAM.panels = [4 1];                              % panels per block
PARAM.rotate = [thetaCone thetaCone thetaCone thetaCone 0];
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) round(rBubble*pi*nPerLenght)];                  % number of element per panel
PARAM.geometryPanel = [0 1 0 1 1];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi];               % theta starting point for the arc
PARAM.rArc = [nan h/2 nan h/2 rBubble];               % theta starting point for the arc
PARAM.x0_Circle = [nan 0 nan L xcmBubble];
PARAM.y0_Circle = [nan r nan r 0];
PARAM.xCrotate = [0 0 0 0 0];
PARAM.yCrotate = [r r r r 0];

%numerics parameters for Laplace
PARAM.typeBClaplace = [2 2 2 2 1];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0 0];     % 0 is straight, 1 is curved (spline)

%numerics parameters for Stokes
PARAM.typeBCstokes = [1 1 1 1 5];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 0 0];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [1 1 1 1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0 1];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [1 1];                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [1 1];
PARAM.deflationConstant = [4*pi*rBubble^2 4*pi*rBubble^2];
PARAM.addFlow = 0;

%repulsive forces
PARAM.repulsiveForces = [1 1];
PARAM.repulsiveOn = 1e-1;
PARAM.coeffRepulsive = 1e5;

%remesh
PARAM.remeshType = [0 0 1 1 1];         % 1 element is split into two parts when coming close to another block, 2 uses a distribution
PARAM.remeshStep = 2;
PARAM.remeshProximity = {[] [] 5 5 3};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = [0 0 L/PARAM.n(3) pi*h/PARAM.n(4) pi*rBubble/PARAM.n(5)];
PARAM.distActivateRemesh = [1 1 1 1 1];
PARAM.adaptCoeff = [4 4 4 4 4];
PARAM.minSizeElemRemesh = [1e-3 1e-3 1e-3 1e-3 1e-3]/5;
PARAM.coeffDist = 2;    % how much smaller than distnce it has to be

%function profile for BCs
PARAM.fluxBC = {0 0 -1 0 []};
PARAM.velBC = {0 0 0 0 []};
PARAM.stressBC{5} = 0;

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)

%build geometry
[initial,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

%create filename
PARAM.filename = chooseFilenameConicalMotor(PARAM,nPerLenght,beta,K,L,Tend,rBubble,xcmBubble,thetaCone,dt,ODE,tol,bubbleCycles);

%remesh function
remeshFun = @(t,x,y,tParametric,PARAM) remeshPanelsFunction(t,x,y,tParametric,PARAM,tParametricBase,nPerLenght);

%event function
event = @(t,x,y,PARAM) eventBlocksCompenetration(t,x,y,PARAM);

%function for motor velocity
fMotor = @(t,x,y,PARAM) computeVelocityConicalMotor(t,x,y,PARAM,K,beta);

tic

%time stepping
[T,X,Y,U,V,yLaplace,yStokes] = RK2panels(fMotor,Tsave,initial,dt,dt,remeshFun,event);

simulationTime = toc;

%save results
display('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

display('The end')



