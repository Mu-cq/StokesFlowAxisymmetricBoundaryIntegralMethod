%compute the motion of a conical motor with diffusion of chemical species

function main_verticalPipeDropParam(alpha,nDrop,nUpload,Bond,lambda,BondUpload,results,uploadRes,Tend,dt,dtUP,ODE,volTol,ODEtol,adaptLevel,remeshTypeBubble,repOn,L,dropFrame,xStart0,TendUP,repUP,ODEup)
    
%results
PARAM.res = results;
PARAM.here = pwd;

%geometrical and physical parameters
%L = 20;                         % channel lenght
H = 1;                          % chennel height
hFilm = 0.2;                    % chennel height
PARAM.visc(4) = lambda;         % viscosity ratio
PARAM.Bond = Bond;              % capillary number
PARAM.D(4) = 0;                 % droplet size
%xStart0 = -L/2+10;

%upload results
BoUp = BondUpload;
OdeUp = ODEup;
x0Up = xStart0;
alphaUP = alpha;
lambdaUP = lambda;
LUP = L;
nUP = nUpload;
%repUP = repOn;
%TendUP = Tend;

%options
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
PARAM.cfunction = 0;
PARAM.STstokes = 1;

%frame of reference
PARAM.dropFrame = dropFrame;    % 0 is lab frame, 1 id drop frame

%time discretization
Tstart = 0;
tol = ODEtol;

%saving options
SaveHowMany = 1e2;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time
PARAM.SaveDataIte = 1;

%geometry parameters
nPerLenght = 10;
PARAM.panels = [3 1];                              % panels per block
PARAM.rotate = [0 0 0 0]/180*pi;
PARAM.n = [round(H*nPerLenght) round(L*nPerLenght) round(H*nPerLenght) nDrop];                  % number of element per panel
PARAM.geometryPanel = [0 0 0 2];                % 0 is a straight line, 1 ia an arc, 2 is a defromable object
PARAM.xStart = [-L/2 -L/2 L/2 nan];             % x starting point for the straight lines
PARAM.xEnd = [-L/2 L/2 L/2 nan];               % x ending point for the straight lines
PARAM.yStart = [0 H H nan];             % y starting point for the straight lines
PARAM.yEnd = [H H 0 nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan nan nan 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan nan nan pi];               % theta starting point for the arc
PARAM.rArc = [nan nan nan alpha];               % theta starting point for the arc
PARAM.x0_Circle = [nan nan nan xStart0];
PARAM.y0_Circle = [nan nan nan 0];
PARAM.xCrotate = [0 0 0 0];
PARAM.yCrotate = [0 0 0 0];
PARAM.ellipseShape = [0 0 0 1];

%numerics parameters for Stokes
PARAM.typeBCstokes = [0 0 0 7];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 1];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 1];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 2];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [0 0 0 2];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [0 2];                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [0 0];
PARAM.deflationConstant = [0 4*pi];
PARAM.addFlow = 0;
PARAM.blockVolCorr = [0 0];
if ODE~=0
    PARAM.blockVolCorr(2) = 1;
end
if PARAM.visc(4)<0.1
    PARAM.Qsource(2) = 0;
    PARAM.deflationBlock(2) = 1;
end

%repulsive forces
PARAM.repulsiveForces = repOn;
PARAM.repulsiveOn = 2e-1;
PARAM.coeffRepulsive = -10;
PARAM.smoothingRep = 1;

%remesh
PARAM.remeshType = [0 1 0 remeshTypeBubble];         % 1 element is split into two parts when coming close to another block, 2 uses a distribution, 4 remesh with a certain distributio (usually for deformable objects)
PARAM.remeshStep = 1;
PARAM.remeshProximity = {4 4 4 2};    % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = 1.1*[H/PARAM.n(1) L/PARAM.n(2) H/PARAM.n(3) alpha*pi/PARAM.n(4)];
PARAM.distActivateRemesh = [1 1 1 1];
PARAM.adaptCoeff = [4 4 4 4];
PARAM.minSizeElemRemesh = [1e-3 1e-3 1e-3 1e-3];
PARAM.coeffDist = 1;    % how much smaller than distnce it has to be
PARAM.distr = [nan nan nan 2];  PARAM.adaptDistr = [nan nan nan adaptLevel];
PARAM.maxNumberTotalElem = 1e3;

%function profile for BCs
PARAM.velBC = {0 0 0};
PARAM.velBCaxial = {0 0 0};
PARAM.velBCradial = {0 0 0};
PARAM.stressBC{4} = 1;

%print to screen
printToScreenStokes(PARAM)
printToScreenTimeStepping(ODE,dt,tol)
printToScreenRemesh(PARAM)
printToScreenVerticalTube(L,PARAM.x0_Circle(4),alpha,PARAM.Bond,PARAM.visc(4),PARAM.dropFrame)

%build geometry
[xInitial,yInitial,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

%create filename
PARAM.filename = chooseFilenameVerticalTube(PARAM,sum(PARAM.n),L,Tend,dt,alpha,xStart0,PARAM.Bond,ODE,PARAM.visc(4),sum(PARAM.repulsiveForces));

%build initial shape
if uploadRes==0
    [xInitial{4},yInitial{4}] = drawBubbleInChannel(xStart0,0,alpha,hFilm,H,nDrop);
elseif uploadRes==1
    [xInitial{4},yInitial{4},PARAM] = uploadBubbleShape(BoUp,OdeUp,x0Up,alphaUP,lambdaUP,LUP,nUP,TendUP,dtUP,sum(repUP),PARAM);
    xcm = center_mass(xInitial{4},yInitial{4});
    xInitial{4} = xInitial{4}+xStart0-xcm;
else
    error('Not implemented')
end

%initial condition
initialXY = zeros(2*numel(xInitial{4}),1);
initialXY(1:2:end-1) = xInitial{4};
initialXY(2:2:end) = yInitial{4};

%function for drop velocity
fRising = @(t,var) computeVelocityRising(t,var,tParametricBase,PARAM);

%remesh function
remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);

%volume correction function
volCorrFunction = @(t,var) volCorrBubble(t,var,4/3*pi*alpha^3);

%output functions and events
if ODE~=0
    outFun = @(t,var,flag) eventBlocksRisingInTube(t,var,flag,tParametricBase,fRising,PARAM);
    if remeshTypeBubble==4
        eventRemesh = @(t,var) activateRemeshVolCorrDistrPanels(t,var,PARAM,4,Tsave,4/3*pi*alpha^3,volTol);
    elseif remeshTypeBubble==5
        eventRemesh = @(t,var) activateRemeshVolCorrProximity(t,var,PARAM,4,Tsave,4/3*pi*alpha^3,volTol);
    else
        error('Not implemented')
    end
    options = odeset('RelTol',tol,'AbsTol',1e-3*tol,'OutputFcn',outFun,'Stats','off','MaxStep',dt,'Events',eventRemesh,'InitialStep',dt*1e-2);
else
    %output, saving and event function
    outFun = @(t,var,T,Y,V) eventBlocksRisingInTubeRK2(t,var,T,Y,V,tParametricBase,PARAM);
end

tic

%time stepping
if ODE==0
    [T,Y,V] = RK2mostGeneralBubbleVolCorr(fRising,Tsave,initialXY,dt,dt,outFun,remeshFunction,volCorrFunction);
elseif ODE~=0
    [T,Y] = odeMatlabRemeshVolCorrBubble(fRising,Tsave,initialXY,outFun,remeshFunction,volCorrFunction,options,ODE,PARAM);
end

simulationTime = toc;

%clear function handle
clear outfun
clear fRising
clear remeshFunction
if ODE~=0
    clear eventRemesh
end

%save results
display('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

display('The end')


