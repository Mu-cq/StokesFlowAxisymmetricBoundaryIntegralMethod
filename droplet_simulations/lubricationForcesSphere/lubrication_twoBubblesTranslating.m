%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%physical parameters
rWall = 1;                      % average cone radius
rBubble1 = 1;                    % radius of the bubble
rBubble2 = 1;                    % radius of the bubble
L = 10;                         % motor lenght
h = 0.2;                        % wall thickness
%manyNdistance = 3*ones(10,1);
%manyDistance = 3:2:50;
manyDistance = 2+logspace(-2,0,10);
rDot1 = -1;
rDot2 = 1;
D1 = 0;
D2 = 0;
%manyN = logspace(1,2,10);
manyN = 200*ones(numel(manyDistance),1);

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

%initialize
errBubble1 = zeros(numel(manyN),1);
errBubble2 = zeros(numel(manyN),1);
F1 = zeros(numel(manyN),1);
F2 = zeros(numel(manyN),1);
F1an = zeros(numel(manyN),1);
F2an = zeros(numel(manyN),1);

for i = 1:numel(manyN)
    
distance = manyDistance(i); 
    
display([num2str(i) ' of ' num2str(numel(manyN))])

%geometry parameters
PARAM.panels = [1 1];                              % panels per block
PARAM.n = [round(manyN(i)) round(manyN(i))];                  % number of element per panel
PARAM.typePanel = [1 1];                % 0 is a straight line, 1 ia an arc

PARAM.kernelFreeSpace = 1;
PARAM.STstokes = 1;
PARAM.addFlow = 0;
PARAM.cfunction = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = [0 0];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity, 4 is prescibed axial velocity , 5 is prescribed normal velocity and tangent stress, 6 is rigid body motion
PARAM.orderVariableStokes = [0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [2 2];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = [0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [1 1];
PARAM.deflationConstant = [4*pi*rBubble1^2 4*pi*rBubble2^2];
PARAM.repulsiveForces = [0 0];
PARAM.rigidBodyMotion = [];           % ridid body motion

%build geometry
[x1,y1] = ellipseCartesian(linspace(0,pi,PARAM.n(1)+1),D1);
x1 = rBubble1*x1 + distance/2;
y1 = rBubble1*y1;
[x2,y2] = ellipseCartesian(linspace(0,pi,PARAM.n(2)+1),D2);
x2 = rBubble2*x2 - distance/2;
y2 = rBubble2*y2;
x{1} = x1;  x{2} = x2;
y{1} = y1;  y{2} = y2;

%plot conical motor shape
figure(1)
plot(x1,y1,'k')
hold on
plot(x2,y2,'k')
plot(x2,-y2,'k')
plot(x1,-y1,'k')
grid on
xlabel('x')
ylabel('r')
axis equal
title('Geometry')
drawnow
hold off

%print to screen
printToScreenStokes(PARAM)

%function profile for BCs
PARAM.velBCaxial{1} = rDot1;
PARAM.velBCradial{1} = 0;
PARAM.velBCaxial{2} = rDot2;
PARAM.velBCradial{2} = 0;
% PARAM.velBC{1} = rDot1;
% PARAM.velBC{2} = rDot2;
% PARAM.stressBC{1} = 0;
% PARAM.stressBC{2} = 0;
% [nx1,ny1] = computeNormalVector(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
% [nx2,ny2] = computeNormalVector(x{2},y{2},PARAM.orderVariableStokes(2),PARAM.orderGeometryStokes(2),PARAM.SPlinesType(2));
% PARAM.velBCaxial{1} = nx1*rDot1;
% PARAM.velBCaxial{2} = nx2*rDot2;
% PARAM.velBCradial{1} = ny1*rDot1;
% PARAM.velBCradial{2} = ny2*rDot2;

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%compute force on bubble
weight1 = integrationOnLineWeightAxis(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
weight2 = integrationOnLineWeightAxis(x{2},y{2},PARAM.orderVariableStokes(2),PARAM.orderGeometryStokes(2),PARAM.SPlinesType(2));
fx1 = yStokes(1:2:2*PARAM.n(1)-1);
fx2 = yStokes(2*PARAM.n(1)+1:2:end-1);
F1(i) = weight1*fx1;
F2(i) = weight2*fx2;

display(['F1=' num2str(F1(i))])
display(['F2=' num2str(F2(i))])

%compute anaytical solution
%cd('analyticalTwoBubbleSebastian')
%[U1an,U2an]=velocity_computation_new(rBubble1,rBubble2,distance,rDot1,rDot2);
%cd ..

%compute force with lubrication
if sum(PARAM.typeBCstokes)==0
    gap = distance/2-rBubble1;
    F1an(i) = -3/2*pi*rDot1*rBubble1^2/gap;
    F2an(i) = -3/2*pi*rDot2*rBubble2^2/gap;
    F1anLog(i) = -3/2*pi*rDot1*rBubble1^2/gap + 27/10*pi*log(1/gap/2);
    F2anLog(i) = -3/2*pi*rDot2*rBubble2^2/gap + 27/10*pi*log(1/gap/2);
elseif sum(PARAM.typeBCstokes)==10
    F1an(i) = 0;
    F2an(i) = 0;
end

display(['Fb1an=' num2str(F1an(i))])
display(['Fb2an=' num2str(F2an(i))])

errBubble1(i) = abs(F1(i)-F1an(i))/abs(F1(i));
errBubble2(i) = abs(F2(i)-F2an(i))/abs(F2(i));

display(['err1=' num2str(errBubble1(i))])
display(['err2=' num2str(errBubble2(i))])

end

if i>1
if manyN(1)~=manyN(2)
    
    figure
    subplot(2,1,1)
    loglog(manyN,errBubble1,'-x')
    grid on
    xlabel('n')
    ylabel('errU1')
    title(['R_1=' num2str(rBubble1) ' R_2=' num2str(rBubble2) ' d=' num2str(distance)])
    
    subplot(2,1,2)
    loglog(manyN,errBubble2,'-x')
    grid on
    xlabel('n')
    ylabel('errU2')
    
elseif manyDistance(1)~=manyDistance(2)
    
figure
    loglog((manyDistance-2*rBubble1)/2,F1,'-x')
    grid on
    xlabel('b')
    ylabel('F')
    hold on
    loglog((manyDistance-2*rBubble1)/2,F1an)
    loglog((manyDistance-2*rBubble1)/2,F1anLog)
    grid on
    legend('F numerical','F analytical','Location','Best')
    title('Force acting on the sphere')
    
    figure
    loglog((manyDistance-2*rBubble1)/2,errBubble1,'-x')
    grid on
    xlabel('b')
    ylabel('errF')
    title('Relative error')
    
end
end

% %plot stresses on wall
% fx = yStokes(1:2:end-3);
% fy = yStokes(2:2:end-2);
% figure
% plot(Xsing,fx)
% hold on
% grid on
% plot(Xsing,fy)
% xlabel('x')
% ylabel('r')
% title('Stress on wall')

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








