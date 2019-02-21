%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%results
PARAM.res = '../results';
PARAM.here = pwd;

%geometrical and physical parameters
PARAM.visc(1) = 0;          % viscosity ratio
alpha = 1;                  % droplet size
Ca = 4*pi*alpha^2*[1];                   % capillary number
PARAM.D(1) = 0;                  % droplet size
nDrop = 400;                 %number of elements
PARAM.ellipseShape = 1;
%manyGap = logspace(-5,-1,10);                    % gap size
manyGap = 0.1;                    % gap size

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
PARAM.geometryPanel = 1;                % 0 is a straight line, 1 ia an arc, 2 is a defromable object
PARAM.xStart = nan;             % x starting point for the straight lines
PARAM.xEnd = nan;               % x ending point for the straight lines
PARAM.yStart = nan;             % y starting point for the straight lines
PARAM.yEnd = nan;              % y ending point for the straight lines
PARAM.thetaStart = 0;               % theta starting point for the arc
PARAM.thetaEnd = pi;               % theta starting point for the arc
PARAM.rArc = alpha;               % theta starting point for the arc
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

%function profile for BCs
PARAM.velBC = {1};
PARAM.stressBC = {1};

%print to screen
printToScreenStokes(PARAM)

%loop on gaps
gapVel = zeros(numel(manyGap),1);
for i = 1:numel(manyGap)
    
    gap = manyGap(i);
    
    printToScreenBubbleCloseToWall(gap,Ca,PARAM.visc)

    %build geometry
    PARAM.x0_Circle = gap+alpha;
    xx2 = linspace(0,1,PARAM.n(1)+1);
    temp = (exp(8*xx2)-1)/(exp(8)-1)*pi;
    temp = flip(diff(temp));
    tParametric{1} = [0 cumsum(temp)];
    [xInitial,yInitial] = buildGeometryPanelsParametric(tParametric,PARAM);
    %[xInitial,yInitial,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);

    %initial condition
    initialXY = zeros(2*numel(xInitial{1}),1);
    initialXY(1:2:end-1) = xInitial{1};
    initialXY(2:2:end) = yInitial{1};

    %function for drop velocity
    Uinterface = computeVelocityBubble(0,initialXY,PARAM);
    
    %gap velocity
    gapVel(i) = Uinterface(end-1);
    
end

if numel(manyGap)>1
    
    figure
    loglog(manyGap,gapVel./manyGap','o-k')
    xlabel('$\epsilon$','interpreter','latex')
    ylabel('$\dot{d}/\epsilon$','interpreter','latex')
    grid on
    
else
    
    disp(['Gap velocity is gapVel=' num2str(gapVel)])

end


%end of simulation
disp('The end')


