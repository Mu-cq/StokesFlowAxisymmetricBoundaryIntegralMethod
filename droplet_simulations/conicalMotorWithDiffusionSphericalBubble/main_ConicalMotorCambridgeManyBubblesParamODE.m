%compute the motion of a conical motor with diffusion of chemical species

function main_ConicalMotorCambridgeManyBubblesParamODE(dt,tol,ODE,thetaCone,BN,bubbleCycles,results)

%results
PARAM.res = results;
PARAM.here = pwd;

%physical parameters
r = 1;                      % radial coordinate for pivoting angle
L = 10;                     % motor lenght
h = 0.2;                    % wall thickness
Hcc = 1.85;
beta = 1.4;
bubbleNucleates = BN;
%bubbleCycles = 15;

%options
%PARAM.newNucleation = 1.5;
PARAM.distCritNucleation = 5;
PARAM.cfunction = 0;
PARAM.STstokes = 1;
PARAM.STlaplace = 1;

%time discretization
Tstart = 0;
Tend = 1;
%dt = 1e-3;
%ODE = 2;
%tol = 1e-3;

%saving options
SaveHowMany = 400;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time

%geometry parameters
nPerLenght = 5;
PARAM.panels = [4 1];                              % panels per block
PARAM.rotate = [thetaCone thetaCone thetaCone thetaCone 0];
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) nan];                  % number of element per panel
PARAM.geometryPanel = [0 1 0 1 1 1];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan nan];             % y starting point for the straight lines
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
PARAM.typeBCstokes = [1 1 1 1 5 5];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 0 0 0];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [1 1 1 1 1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0 1 1];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [1 1];                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [1 1 1];
PARAM.deflationConstant = [4*pi 4*pi 4*pi];
PARAM.addFlow = 0;

%repulsive forces
PARAM.repulsiveForces = [1 1 1];
PARAM.repulsiveOn = 1e-1;
PARAM.debeyLenght = 45;
PARAM.coeffRepulsive = 1e6*sin(thetaCone);

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
PARAM.filename = chooseFilenameConicalMotor(PARAM,nPerLenght,beta,Hcc,bubbleNucleates,L,Tend,thetaCone,dt,ODE,tol,bubbleCycles);

%initial condition
initial = [0 PARAM.x0_Circle(5) PARAM.rArc(5)]';
        
tic

%start bubble cycles
manyT = cell(bubbleCycles,1);
manyY = cell(bubbleCycles,1);
for i = 1:bubbleCycles
    
    display(['Bubble cycle number ' num2str(i)])
    
    %number of bubbles
    nBubble = numel(PARAM.n)-4;
    
    %function for motor velocity, one bubble
    f1bubble = @(t,var) computeVelocityConicalMotorODE(t,var,tParametricBase,nPerLenght,PARAM,Hcc,beta);

    %function for motor velocity, two bubble
    f2bubbles = @(t,var) computeVelocityConicalMotorTwoBubblesODE(t,var,tParametricBase,nPerLenght,PARAM,Hcc,beta);
    
    if nBubble==1
        
        %function to solve
        fToSolve = f1bubble;

        %output and event function
        outFun = @(t,var,flag) eventBlocksODE(t,var,flag,tParametricBase,f1bubble,PARAM);
    
    elseif nBubble==2
        
        %function to solve
        fToSolve = f2bubbles;
        
        %output and event function
        outFun = @(t,var,flag) eventBlocksTwoBubblesODE(t,var,flag,tParametricBase,f2bubbles,beta,Hcc,PARAM);
        
    else
        
       error('Current implementation supports maximum two bubbles')
        
    end
    
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
    [breakNow,initial,tParametricBase,Tsave,PARAM] = initialConditionNextCycle(tParametricBase,nBubble,manyT{i},manyY{i},fToSolve,Tsave,Tend,SaveHowMany,beta,Hcc,PARAM);
    
    if breakNow==1
        
       break; 
        
    end

end

simulationTime = toc;

%save results
display('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

display('The end')


