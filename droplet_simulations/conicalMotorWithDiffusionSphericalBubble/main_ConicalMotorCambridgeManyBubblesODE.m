%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/droplet_simulations/results';
PARAM.here = pwd;

warning('Careful with repulsive forces')

%physical parameters
r = 1;                          % radial coordinate for pivoting angle
L = 10;                         % motor lenght
h = 0.2;                        % wall thickness
thetaCone = 2/180*pi;           % inclinitation of the cone
Hcc = 0.1;
beta = 1;
bubbleNucleates = 20;
bubbleCycles = 100;
maxNbubbles = 10;

%options nucleation
%PARAM.newNucleation = 1.5;
endCycle = [];   %if 1 stops when center of mass exit the cone, if 2 stops when second bubble nucleates, if
PARAM.maxRadiusNucleation = 1;
PARAM.smallestRadius = 0;
PARAM.distCritNucleation = 5;

%options numerics
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.STlaplace = 1;

%time discretization
Tstart = 0;
Tend = 1;
dt = 1e-3;
ODE = 2;
tol = 1e-3;

%saving options
SaveHowMany = 100;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time

%geometry parameters
nPerLenght = 5;
PARAM.panels = [4 1];                              % panels per block
PARAM.rotate = [thetaCone thetaCone thetaCone thetaCone 0];
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) nan];                  % number of element per panel
PARAM.geometryPanel = [0 1 0 1 ones(1,maxNbubbles)];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0*ones(1,maxNbubbles)];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi*ones(1,maxNbubbles)];               % theta starting point for the arc
PARAM.rArc = [nan h/2 nan h/2];               % theta starting point for the arc
PARAM.x0_Circle = [nan 0 nan L nan];
PARAM.y0_Circle = [nan r nan r 0];
PARAM.xCrotate = [0 0 0 0 0];
PARAM.yCrotate = [r r r r 0];

%numerics parameters for Laplace
PARAM.typeBClaplace = [2 2 2 2 ones(1,maxNbubbles)];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0 0*ones(1,maxNbubbles)];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0 0*ones(1,maxNbubbles)];     % 0 is straight, 1 is curved (spline)

%numerics parameters for Stokes
PARAM.typeBCstokes = [1 1 1 1 5*ones(1,maxNbubbles)];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0 0*0*ones(1,maxNbubbles)];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0 0*0*ones(1,maxNbubbles)];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 0 0*0*ones(1,maxNbubbles)];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [1 1 1 1 ones(1,maxNbubbles)];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0 1*ones(1,maxNbubbles)];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [1 1];                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = ones(1,maxNbubbles+1);
PARAM.deflationConstant = 4*pi*ones(1,maxNbubbles+1);
PARAM.addFlow = 0;

%repulsive forces
PARAM.repulsiveForces = ones(1,10);
PARAM.repulsiveOn = 1e-1;
PARAM.debeyLenght = 45;
PARAM.coeffRepulsive = 1e7*sin(thetaCone);

%remesh
PARAM.remeshType = [0 0 1 1 ones(1,maxNbubbles)];         % 1 element is split into two parts when coming close to another block, 2 uses a distribution
PARAM.remeshStep = 2;
PARAM.remeshProximity = {[] [] 5 5 3 3};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = [0 0 L/PARAM.n(3) pi*h/PARAM.n(4) nan*ones(1,maxNbubbles)];
PARAM.distActivateRemesh = [1 1 1 1 ones(1,maxNbubbles)];
PARAM.adaptCoeff = [4 4 4 4 4*ones(1,maxNbubbles)];
PARAM.minSizeElemRemesh = [1e-3 1e-3 1e-3 1e-3 1e-3*ones(1,maxNbubbles)]/2;
PARAM.coeffDist = 2;    % how much smaller than distnce it has to be

%function profile for BCs
PARAM.fluxBC = {0 0 -1 0 []};
PARAM.velBC = {0 0 0 0 []};
PARAM.stressBC{5} = 0;
PARAM.stressBC{6} = 0;

%choose bubble position and size
PARAM = chooseBubblePositionAndSize(PARAM,beta,Hcc,nPerLenght,bubbleNucleates);

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)
printToScreenTimeStepping(ODE,dt,tol)

%build geometry
[~,~,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

%create filename
PARAM.filename = chooseFilenameConicalMotor(PARAM,nPerLenght,beta,Hcc,bubbleNucleates,L,Tend,thetaCone,dt,ODE,tol,bubbleCycles,endCycle);

%initial condition
initial = [0 PARAM.x0_Circle(5) PARAM.rArc(5)]';
        
tic

%start bubble cycles
manyT = cell(bubbleCycles,1);
manyY = cell(bubbleCycles,1);
for i = 1:bubbleCycles
    
    disp(['Bubble cycle number ' num2str(i) ' begins with ' num2str(numel(PARAM.n)-4) ' bubbles'])
    
    if numel(PARAM.n)-4>maxNbubbles
        disp('Too many bubbles')
    end
    
    %function for motor velocity, two bubble
    fToSolve = @(t,var) computeVelocityConicalMotorManyBubblesODE(t,var,tParametricBase,nPerLenght,PARAM,Hcc,beta);
    
    %output and event function
    outFun = @(t,var,flag) eventBlocksManyBubblesODE(t,var,flag,tParametricBase,fToSolve,beta,Hcc,PARAM);
    
    %ODE options
    options = odeset('RelTol',tol,'AbsTol',1e-6,'OutputFcn',outFun,'Stats','off','MaxStep',dt);

    %time stepping
    if ODE==0
        [manyT{i},manyY{i}] = RK2mostGeneral(fToSolve,Tsave,initial,dt,dt,outFun);
    elseif ODE==1
        [manyT{i},manyY{i}] = ode45(fToSolve,Tsave, initial, options);
    elseif ODE==2
        [manyT{i},manyY{i}] = ode23t(fToSolve,Tsave, initial, options);
    elseif ODE==3
        [manyT{i},manyY{i}] = ode23s(fToSolve,Tsave, initial, options);
    end
    
    %prepare initial condition next cycle
    [breakNow,initial,tParametricBase,Tsave,PARAM] = initialConditionNextCycleManyBubbles(tParametricBase,manyT{i},manyY{i},Tsave,Tend,SaveHowMany,outFun,PARAM);
    
    if breakNow==1
        
       break;
        
    end

end

simulationTime = toc;

%save results
disp('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

disp('The end')


