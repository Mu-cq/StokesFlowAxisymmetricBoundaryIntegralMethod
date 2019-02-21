%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
rBubble0 = 1;                    % radius of the bubble
gap0 = 1;
%manyGap = 1e-3;
rDot = 1;
D = 0;
elem = 200;

%time stepping
dt = 0.01;
loop = 5e3;
T = 0:dt:loop*dt;

%options
plotShape = 0;

%geometry parameters
PARAM.panels = 1;                              % panels per block
PARAM.n = round(elem);                  % number of element per panel
PARAM.geometryPanel = 1;                % 0 is a straight line, 1 ia an arc
PARAM.xStart = nan;             % x starting point for the straight lines
PARAM.xEnd = nan;               % x ending point for the straight lines
PARAM.yStart = nan;             % y starting point for the straight lines
PARAM.yEnd = nan;               % y ending point for the straight lines
PARAM.thetaStart = 0;               % theta starting point for the arc
PARAM.thetaEnd = pi;               % theta starting point for the arc
PARAM.rArc = rBubble0;               % theta starting point for the arc
PARAM.y0_Circle = 0;

%options
PARAM.STstokes = 1;
PARAM.kernelFreeSpace = 2;  PARAM.posWall = 0;
PARAM.addFlow = 0;
PARAM.cfunction = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = 0;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 0;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 0;    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = 2;            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = 1;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = 1;
PARAM.deflationConstant = 4*pi*rBubble0^2;
PARAM.repulsiveForces = 0;

%print to screen
printToScreenStokes(PARAM)

%euler loop
sphereVel = zeros(loop,1);
gap = zeros(loop,1);
rBubble = zeros(loop,1);
gap(1) = gap0;
rBubble(1) = rBubble0;
for i = 1:loop
    
    display([num2str(i) ' of ' num2str(loop)])
    
    %build geometry
    PARAM.x0_Circle = gap(i)+rBubble(i);
    PARAM.rArc = rBubble(i);
    
    %sphere velocity
    [Usphere,exitFunction] = computeVelSphere(plotShape,rDot,PARAM);
    
    %update sphere position and radius
    gap(i+1) = gap(i)+dt*(-rDot+Usphere);
    rBubble(i+1) = rBubble(i)+dt*rDot;
    
    if exitFunction==1
        break;
    end

end

if PARAM.typeBCstokes==0
    typeBC = 'rigid shell';
elseif PARAM.typeBCstokes==5
    typeBC = 'bubble';
end

figure
semilogy(T(1:i+1),gap(1:i+1),'k')
[minGap,indMin] = min(gap(1:i+1));
hold on
semilogy(T(indMin),minGap,'.r','MarkerSize',30)
xyLabelTex('t','d')
grid on
title(typeBC,'interpreter','latex')

figure
semilogy(rBubble0+rDot*T(1:i+1),gap(1:i+1),'k')
[minGap,indMin] = min(gap);
hold on
plot(rBubble0+rDot*T(indMin),minGap,'.r','MarkerSize',30)
xyLabelTex('a','d')
grid on
title(typeBC,'interpreter','latex')

function [Usphere,status] = computeVelSphere(plotShape,rDot,PARAM)

    status = 0;

    xx2 = linspace(0,1,PARAM.n(1)+1);
    cluster = 16;
    temp = (exp(cluster*xx2)-1)/(exp(cluster)-1)*pi;
    temp = flip(diff(temp));
    tParametric{1} = [0 cumsum(temp)];
    [x,y,PARAM.minSize,PARAM.maxSize] = buildGeometryPanelsParametric(tParametric,PARAM);
    
    %function profile for BCs
    temp = tParametric{1};  temp = (temp(1:end-1)+temp(2:end))/2;
    if PARAM.typeBCstokes(1)==0
        PARAM.velBCaxial{1} = rDot*cos(temp);
        PARAM.velBCradial{1} = rDot*sin(temp);
    elseif PARAM.typeBCstokes(1)==1
        PARAM.velBC{1} = rDot;
    elseif PARAM.typeBCstokes(1)==5
        PARAM.velBC{1} = rDot;
        PARAM.stressBC{1} = 0;
    end
    
    if plotShape==1
        figure(1)
        plot(x{1},y{1},'xk')
        hold on
        plot(x{1},-y{1},'k')
        grid on
        xlabel('x')
        ylabel('r')
        axis equal
        if PARAM.typeBCstokes(1)==5
            title('Bubble')
        else
            title('Solid sphere')
        end
        drawnow
        hold off
    end
    
    if min(x{1})<0
        warning('Bubble interface crosses the wall')
        status = 1;
    end

    %solve Stokes equation
    yStokes = BEM_Stokes(x,y,PARAM);
    
    Usphere = yStokes(end);

end






