%compute the motion of a conical motor with diffusion of chemical species

function main_oneDropBEMparam(nDrop,Ca,visc,dt,Tend,curved,res,nameShort)

%results
PARAM.res = res;
PARAM.here = pwd;

%choose algorithm
PARAM.algorithm = 1;    %1 is time marching, 2 is minimal seed

%geometrical and physical parameters
PARAM.visc(1) = visc;                  % viscosity ratio
PARAM.Ca = Ca;                     % capillary number
PARAM.Bond = [];                    % bond number
xStart0 = 0;
V0 = 4/3*pi;
energyPerturb = [];

%initial shape
PARAM.ellipseShape = 1;
PARAM.D(1) = 0;                  % droplet deformation
PARAM.f{1} = [-0.5 0 1];    PARAM.whichF{1} = [2 4 10];

%upload results
PARAM.uploadRes = 0;
PARAM.BoUp = [];
PARAM.CaUP = PARAM.Ca;
PARAM.OdeUp = 0;
PARAM.x0Up = 0;
PARAM.lambdaUP = PARAM.visc(1);
PARAM.nUP = 10;
PARAM.TendUP = 100;
PARAM.dtUP = 0.0;

%options
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];

%frame of reference
PARAM.dropFrame = 1;    % 0 is lab frame, 1 id drop frame

%time discretization
Tstart = 0;
%Tend = 100;
%dt = 1e-1;
PARAM.ODE = 0;
tol = 1e-3;

%saving options
SaveHowMany = 1e2;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time
PARAM.SaveDataIte = 0;

%geometry parameters
%nDrop = 20;
PARAM.panels = 1;                              % panels per block
PARAM.rotate = [0 0 0 0]/180*pi;
PARAM.n = nDrop;                  % number of element per panel
PARAM.geometryPanel = 2;                % 0 is a straight line, 1 ia an arc, 2 is a defromable object
PARAM.xStart = nan;             % x starting point for the straight lines
PARAM.xEnd = nan;               % x ending point for the straight lines
PARAM.yStart = nan;             % y starting point for the straight lines
PARAM.yEnd = nan;               % y ending point for the straight lines
PARAM.thetaStart = nan;               % theta starting point for the arc
PARAM.thetaEnd = nan;               % theta starting point for the arc
PARAM.rArc = nthroot(3/4/pi*V0,3);               % theta starting point for the arc
PARAM.x0_Circle = xStart0;
PARAM.y0_Circle = 0;
PARAM.xCrotate = 0;
PARAM.yCrotate = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = 2;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 1;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = curved;    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = 2;            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = 2;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = 0;           % 0 is not inflating, 1 is inflating
PARAM.blockType = 2;                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = 0;
PARAM.deflationConstant = 0;
PARAM.addFlow = 2;

if PARAM.visc(1)<0.1
    PARAM.deflationBlock = 1;
    PARAM.deflationConstant = [];
    PARAM.Qsource = 0;
end

%repulsive forces
PARAM.repulsiveForces = 0;
PARAM.repulsiveOn = 2e-1;
PARAM.coeffRepulsive = 1e-2;
PARAM.smoothingRep = [];

%remesh
PARAM.remeshType = 4;         % 1 element is split into two parts when coming close to another block, 2 uses a distribution, 4 remesh with a certain distributio (usually for deformable objects)
PARAM.remeshStep = 1;
PARAM.remeshProximity = {};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = 1.1*pi/PARAM.n(1);
PARAM.distActivateRemesh = 1;
PARAM.adaptCoeff = 4;
PARAM.minSizeElemRemesh = 1e-3/2;
PARAM.coeffDist = 1;    % how much smaller than distnce it has to be
PARAM.distr = 0;  PARAM.adaptDistr = 0;
PARAM.maxNumberTotalElem = 1e3;

%function profile for BCs
PARAM.velBC = {};
PARAM.velBCaxial = {};
PARAM.velBCradial = {};
PARAM.stressBC{1} = 1;

%print to screen
printToScreenStokes(PARAM)
printToScreenTimeStepping(PARAM.ODE,dt,tol)
printToScreenRemesh(PARAM)
printToScreenOneDrop(PARAM.Ca,PARAM.Bond,PARAM.visc(1),PARAM.dropFrame)

%create filename
PARAM.filename = chooseFilenameOneDropBEM(PARAM,dt,Tend,[],PARAM.ODE,nameShort);

tic

%run simulation
if PARAM.algorithm==1
    %time stepping
    allRes = runTimeSteppingBEM(Tsave,dt,dt,PARAM,V0);
elseif PARAM.algorithm==2
    %minimal seed
    %allRes = runTimeSteppingBEM(Tsave,dt,dt,PARAM,V0);
    minimalSeedOneDropBEM(Tsave,dt,dt,PARAM,V0,energyPerturb)
end

simulationTime = toc;

%save results
if PARAM.algorithm==1
    disp('Save results')
    cd(PARAM.res)
    save(PARAM.filename)
    cd(PARAM.here)
end

disp('The end')


