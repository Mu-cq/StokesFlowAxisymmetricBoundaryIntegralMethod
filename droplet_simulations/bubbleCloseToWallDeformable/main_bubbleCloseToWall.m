%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%results
PARAM.res = '../results';
PARAM.here = pwd;

%geometrical and physical parameters
gap = 0.1;                    % chennel height
PARAM.visc(1) = 0;          % viscosity ratio
Ca = 100;                   % capillary number
alpha = 1;                  % droplet size
PARAM.D(1) = 0;                  % droplet size
nDrop = 40;                 %number of elements
PARAM.ellipseShape = 1;

%options
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.kernelFreeSpace = 2;  PARAM.posWall = 0;

%time discretization
Tstart = 0;
Tend = 100;
dt = 0.01;
ODE = 0;
tol = 1e-3;
PARAM.Unormal = 0;

%saving options
SaveHowMany = 1e2;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time
PARAM.SaveDataIte = 0;

%geometry parameters
PARAM.panels = 1;                              % panels per block
PARAM.rotate = [0 0 0 0]/180*pi;
PARAM.n = nDrop;                  % number of element per panel
PARAM.geometryPanel = 2;                % 0 is a straight line, 1 ia an arc, 2 is a defromable object
PARAM.xStart = nan;             % x starting point for the straight lines
PARAM.xEnd = nan;               % x ending point for the straight lines
PARAM.yStart = nan;             % y starting point for the straight lines
PARAM.yEnd = nan;              % y ending point for the straight lines
PARAM.thetaStart = 0;               % theta starting point for the arc
PARAM.thetaEnd = pi;               % theta starting point for the arc
PARAM.rArc = alpha;               % theta starting point for the arc
PARAM.x0_Circle = alpha+gap;
PARAM.y0_Circle = 0;
PARAM.xCrotate = 0;
PARAM.yCrotate = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = 2;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 1;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 1;    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = 2;            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = 2;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = 2;           % 0 is not inflating, 1 is inflating
PARAM.blockType = 2;                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = 1;
PARAM.deflationConstant = 0;
PARAM.addFlow = 0;
PARAM.Qsource(1) = Ca;
%PARAM.blockVolCorr = 0;
volTol = 1e-3;

%repulsive forces
PARAM.repulsiveForces = [0 0];
PARAM.repulsiveOn = 2e-1;
PARAM.coeffRepulsive = 1e-2;
PARAM.smoothingRep = [];

%remesh
PARAM.remeshType = 4;         % 1 element is split into two parts when coming close to another block, 2 uses a distribution, 4 remesh with a certain distributio (usually for deformable objects)
PARAM.remeshStep = 1;
PARAM.remeshProximity = {};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = 1.1*alpha*pi/PARAM.n(1);
PARAM.distActivateRemesh = 1;
PARAM.adaptCoeff = 4;
PARAM.minSizeElemRemesh = 1e-3/2;
PARAM.coeffDist = 1;    % how much smaller than distnce it has to be
PARAM.distr = 3;  PARAM.adaptDistr = .5;
PARAM.maxNumberTotalElem = 1e3;

%function profile for BCs
PARAM.velBC = {1};
PARAM.stressBC = {1};

%print to screen
printToScreenStokes(PARAM)
printToScreenTimeStepping(ODE,dt,tol)
printToScreenRemesh(PARAM)
printToScreenBubbleCloseToWall(gap,Ca,PARAM.visc)

%build geometry
[xInitial,yInitial,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

%create filename
PARAM.filename = chooseFilenameBubbleCloseToWall(PARAM,nDrop,Tend,dt,alpha,Ca,gap,ODE,PARAM.visc);

%initial condition
initialXY = zeros(2*numel(xInitial{1}),1);
initialXY(1:2:end-1) = xInitial{1};
initialXY(2:2:end) = yInitial{1};

%function for drop velocity
fBubble = @(t,var) computeVelocityBubble(t,var,PARAM);

%remesh function
remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);
    
%output functions and events
if ODE~=0
    outFun = @(t,var,flag) eventBlocksBubbleCloseToWall(t,var,flag,PARAM);
    eventRemeshVolCorr = @(t,var) activateRemeshVolCorrProximity(t,var,PARAM,4,Tsave,4/3*pi*alpha^3,volTol);
    options = odeset('RelTol',tol,'AbsTol',1e-3*tol,'OutputFcn',outFun,'Stats','off','MaxStep',dt,'Events',eventRemeshVolCorr,'InitialStep',dt*1e-2);
else
    %output, saving and event function
    outFun = @(t,var,T,Y,V) eventBlocksBubbleCloseToWallRK2(t,var,T,Y,V,PARAM);
end

tic

%time stepping
if ODE==0
    [T,Y,V] = RK2mostGeneralBubble(fBubble,Tsave,initialXY,dt,dt,outFun,remeshFunction);
elseif ODE~=0
    [T,Y] = odeMatlabRemeshVolCorrBubble(fBubble,Tsave,initialXY,outFun,remeshFunction,options,ODE,PARAM);
end

simulationTime = toc;

%clear function handle
clear outfun
clear fBubble
clear remeshFunction
clear volCorrFunction
if ODE~=0
    clear eventRemesh
end

%save results
disp('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

disp('The end')


