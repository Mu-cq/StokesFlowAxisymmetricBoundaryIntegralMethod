%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
rWall = 1;                      % average cone radius
rBubble1 = 1;                    % radius of the bubble
rBubble2 = 1;                    % radius of the bubble
manyGap = logspace(-5,-1,10);
rDot1 = 1;
rDot2 = 1;
D1 = 0;
D2 = 0;
elem = 200;

%options
plotShape = 1;

%grid for velocity field
oneline = 2;
if oneline==0
    
    xxx = linspace(-4,4,80);
    yyy = linspace(0,4,40);
    [X,Y] = meshgrid(xxx,yyy);
    
elseif oneline==1
    
    rrr = logspace(0,4,100);
    theta = pi/4;
    X = rrr.*cos(theta);
    Y = rrr.*sin(theta);
    
end

%geometry parameters
PARAM.panels = [1 1];                              % panels per block
PARAM.n = [round(elem) round(elem)];                  % number of element per panel
PARAM.geometryPanel = [1 1];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [nan nan];             % x starting point for the straight lines
PARAM.xEnd = [nan nan];               % x ending point for the straight lines
PARAM.yStart = [nan nan];             % y starting point for the straight lines
PARAM.yEnd = [nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [0 0];               % theta starting point for the arc
PARAM.thetaEnd = [pi pi];               % theta starting point for the arc
PARAM.rArc = [rBubble1 rBubble2];               % theta starting point for the arc
PARAM.y0_Circle = [0 0];

%options
PARAM.STstokes = 1;
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
PARAM.addFlow = 0;
PARAM.cfunction = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = [1 1];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [2 2];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = [1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [1 1];
PARAM.deflationConstant = [4*pi*rBubble1^2 4*pi*rBubble2^2];
PARAM.repulsiveForces = [0 0];

%print to screen
printToScreenStokes(PARAM)

%function profile for BCs
PARAM.velBC{1} = rDot1;
PARAM.stressBC{1} = 0;
PARAM.velBC{2} = rDot2;
PARAM.stressBC{2} = 0;

%euler loop
U1 = zeros(numel(manyGap),1);
U2 = zeros(numel(manyGap),1);
for i = 1:numel(manyGap)
    
    display([num2str(i) ' of ' num2str(numel(manyGap))])
    
    %build geometry
    gap = manyGap(i);
    PARAM.x0_Circle = [gap+rBubble1 -gap-rBubble2];
    xx2 = linspace(0,1,PARAM.n(2)+1);
    tParametric{2} = (exp(8*xx2)-1)/(exp(8)-1)*pi;
    tInv = flip(diff(tParametric{2}));
    tParametric{1} = [0 cumsum(tInv)];
    [x,y,PARAM.minSize,PARAM.maxSize] = buildGeometryPanelsParametric(tParametric,PARAM);

    if plotShape==1
        figure(1)
        plot(x{1},y{1},'xk')
        hold on
        plot(x{2},y{2},'xk')
        plot(x{2},-y{2},'k')
        plot(x{1},-y{1},'k')
        grid on
        xlabel('x')
        ylabel('r')
        axis equal
        if PARAM.typeBCstokes(1)==5
            title('Clean bubble')
        else
            title('Dirty bubble')
        end
        drawnow
        hold off
    end

    %solve Stokes equation
    [yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);
    
    %bubble velocity
    U1(i) = yStokes(end-1);
    U2(i) = yStokes(end);

    %compute force on bubble
    weight1 = integrationOnLineWeightAxis(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
    weight2 = integrationOnLineWeightAxis(x{2},y{2},PARAM.orderVariableStokes(2),PARAM.orderGeometryStokes(2),PARAM.SPlinesType(2));
    if PARAM.typeBCstokes(1)==5
        [nx1,ny1] = computeNormalVector(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
        [nx2,ny2] = computeNormalVector(x{2},y{2},PARAM.orderVariableStokes(2),PARAM.orderGeometryStokes(2),PARAM.SPlinesType(2));
        fx1 = yStokes(1:2:2*PARAM.n(1)-1).*nx1';
        fx2 = yStokes(2*PARAM.n(1)+1:2:end-3).*nx2';
    else
        fx1 = yStokes(1:2:2*PARAM.n(1)-1);
        fx2 = yStokes(2*PARAM.n(1)+1:2:end-3);
    end
    F1 = weight1*fx1;
    F2 = weight2*fx2;

    if abs(F1)>1e-10||abs(F2)>1e-10
        error('Bubble are not force free')
    end

    manyGap(i) = (min(x{1})-max(x{2}))/2;

end

%plot gap
figure
loglog(manyGap,U1-rDot1,'kx--')
xlabel('gap')
ylabel('Ugap')
grid on
title('Velocity scaling')

%compute velocity field
if oneline==1||oneline==0

    %compute velocity field
    [XsingPlot,YsingPlot,ux,uy,p] = computeVelPressField(X,Y,x,y,yStokes,0,PARAM,0,0);

end

if oneline==0

    figure
    streamslice(XsingPlot,YsingPlot,ux,uy)
    hold on
    contourf(XsingPlot,-YsingPlot,sqrt(ux.^2+uy.^2))
    colorbar
    plot(x1,y1,'k')
    plot(x1,-y1,'k')
    plot(x2,y2,'k')
    plot(x2,-y2,'k')
    xlabel('x')
    ylabel('r')
    axis([min(min(XsingPlot)) max(max(XsingPlot)) -max(max(YsingPlot)) max(max(YsingPlot))])
    axis equal
    title('Velocity field')

elseif oneline==1
    
    figure
    loglog(sqrt(XsingPlot.^2+YsingPlot.^2),sqrt(ux.^2+uy.^2),'-x')
    xlabel('d')
    ylabel('|U|')
    grid on
    title('Velocity decay')
    
end








