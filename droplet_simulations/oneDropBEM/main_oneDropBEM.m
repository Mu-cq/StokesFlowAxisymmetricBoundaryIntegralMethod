%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%results
PARAM.res = '../results';
nameShort = '';
PARAM.here = pwd;

%% CHOOSE ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.algorithm = 3;    %1 is time marching, 2 is minimal seed 3 is edge tracking


%% PHYSICS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.visc(1) = 1;                 % viscosity ratio
PARAM.Ca = [];                     % capillary number
PARAM.Bond = 1.5;                    % bond number
xStart0 = 0;
V0 = 4/3*pi;
energyPerturb = [];


%% SHAPE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.ellipseShape = 1;
PARAM.D(1) = 0.875;                  % droplet deformation
PARAM.f{1} = [1 0 0];    PARAM.whichF{1} = [2 4 10];


%% UPLOAD OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.uploadRes = 0;
PARAM.BoUp = [];
PARAM.CaUP = PARAM.Ca;
PARAM.OdeUp = 0;
PARAM.x0Up = 0;
PARAM.lambdaUP = PARAM.visc(1);
PARAM.nUP = 10;
PARAM.TendUP = 100;
PARAM.dtUP = 0.1;


%% NUMERICS OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.cfunction = 1;
PARAM.STstokes = 1;
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];


%% FRAME OF REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.dropFrame = 1;    % 0 is lab frame, 1 id drop frame


%% TIME DISCRETIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tstart = 0;
Tend = 60;
dt = 1e-3;
PARAM.ODE = 0;
PARAM.tol = 1e-3;
initialDT = dt;


%% EDGE TRACKING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.blockEdge = 1;
PARAM.edgeLoop = 1001;
PARAM.deltaEdge = 1e-6;
PARAM.edgeStartingLoop = 19;
PARAM.Dedge = [0.8 0.875];
PARAM.normEdge = 2;
PARAM.bisection = 1;
PARAM.interpOrExtrap = 2;
PARAM.forceIntegrationAfterDT = 5;


%% CONVERGENCE CRITERION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.resConverge = 1e-3;
PARAM.convergeShapeEdge = 0.02;


%% SAVING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SaveHowMany = Tend*20;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time
PARAM.SaveDataIte = 0;
PARAM.saveData = 1;


%% GEOMETRY PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDrop = 70;
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


%% NUMERICS PARAMETERS FOR STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.typeBCstokes = 7;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 1;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 1;    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = 2;            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = 2;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = 0;           % 0 is not inflating, 1 is inflating
PARAM.blockType = 2;                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = 0;
PARAM.deflationConstant = 0;
PARAM.addFlow = 0;
PARAM.Qsource = 0;
%PARAM.blockVolCorr = 0;
volTol = 1e-3;


%% REPULSIVE FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.repulsiveForces = 0;
PARAM.repulsiveOn = 2e-1;
PARAM.coeffRepulsive = 1e-2;
PARAM.smoothingRep = [];


%% REMESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.normRemesh = 1e-4;
PARAM.remeshType = 4;         % 1 element is split into two parts when coming close to another block, 2 uses a distribution, 4 remesh with a certain distributio (usually for deformable objects)
PARAM.remeshStep = 1;
PARAM.remeshProximity = {};    % in case remesh by proximity, indicate which with panel to chek proximity
[~,~,Lside] = ellipseCartesian(linspace(0,pi,100),PARAM.D);
PARAM.maxElem = 1.1*Lside/PARAM.n(1);
PARAM.distActivateRemesh = 1;
PARAM.adaptCoeff = 4;
PARAM.minSizeElemRemesh = 1e-3;
PARAM.coeffDist = 1;    % how much smaller than distnce it has to be
PARAM.distr = 1;  PARAM.adaptDistr = 0.1;
PARAM.maxNumberTotalElem = 1e3;


%% FUNCTION FOR BCS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.velBC = {};
PARAM.velBCaxial = {};
PARAM.velBCradial = {};
PARAM.stressBC{1} = 1;


%% PRINT TO SCREEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printToScreenStokes(PARAM)
printToScreenTimeStepping(PARAM.ODE,dt,PARAM.tol)
printToScreenRemesh(PARAM)


%% FILENAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.filename = chooseFilenameOneDropBEM(PARAM,dt,Tend,[],PARAM.ODE,nameShort);

tic

%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PARAM.algorithm==1
    %time stepping
    allRes = runTimeSteppingBEM(Tsave,dt,initialDT,PARAM,V0);
elseif PARAM.algorithm==2
    %minimal seed
    minimalSeedOneDropBEM(Tsave,dt,dt,PARAM,V0,energyPerturb)
elseif PARAM.algorithm==3
    %edge tracking
    [Tedge,Yedge,Vedge] = edgeTrackingBEM(Tsave,dt,initialDT,PARAM,V0);
end

simulationTime = toc;

%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PARAM.algorithm==1 || PARAM.algorithm==3
    disp('Save results')
    cd(PARAM.res)
    save(PARAM.filename)
    cd(PARAM.here)
end

disp('The end')


