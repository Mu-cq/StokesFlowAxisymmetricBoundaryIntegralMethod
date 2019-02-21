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
L = 10;                         % motor lenght
h = 0.2;                        % wall thickness
rDot = 1;
D1 = 0;
D2 = 0;
%manyN = logspace(1,2,10);
manyN = 200*ones(1,1);

%grid for velocity field
oneline = 0;
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
errSphere = zeros(numel(manyN),1);

for i = 1:numel(manyN)
    
display([num2str(i) ' of ' num2str(numel(manyN))])

%geometry parameters
PARAM.panels = 1;                              % panels per block
PARAM.n = round(manyN(i));                  % number of element per panel
PARAM.typePanel = 1;                % 0 is a straight line, 1 ia an arc

PARAM.STstokes = 1;
PARAM.addFlow = 0;
PARAM.cfunction = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = 0;           % 0 is prescribed axial velocity and radial velocity, 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 0;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 0;    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = 2;            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = 0;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 0;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = 1;
PARAM.deflationConstant = 4*pi*rBubble1^2;
PARAM.repulsiveForces = 0;

%build geometry
[x1,y1] = ellipseCartesian(linspace(0,pi,PARAM.n(1)+1),D1);
x1 = rBubble1*x1;
y1 = rBubble1*y1;
x{1} = x1;
y{1} = y1;

%plot conical motor shape
figure(1)
plot(x1,y1,'k')
hold on
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

%compute vel for BC
[nx,ny] = computeNormalVector(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
axial = nx*rDot;
radial = ny*rDot;

%function profile for BCs
PARAM.velBCaxial{1} = axial;
PARAM.velBCradial{1} = radial;

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%compute force on bubble
weight1 = integationOnLineWeightAxis(x{1},y{1},Ysing(1:PARAM.n(1))',PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
fx = yStokes(1:2:2*PARAM.n(1)-1);
F = weight1*fx;

display(['F=' num2str(F)])

%compute anaytical solution
%cd('analyticalTwoBubbleSebastian')
%[U1an,U2an]=velocity_computation_new(rBubble1,rBubble2,distance,rDot1,rDot2);
%cd ..

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
    
%     figure
%     subplot(2,1,1)
%     plot(manyDistance,errBubble1,'-x')
%     grid on
%     xlabel('d')
%     ylabel('errU1')
%     title(['R_1=' num2str(rBubble1) ' R_2=' num2str(rBubble2) ' n=' num2str(manyN(1))])
%     
%     subplot(2,1,2)
%     plot(manyDistance,errBubble2,'-x')
%     grid on
%     xlabel('d')
%     ylabel('errU2')
    
    figure
    subplot(2,1,1)
    loglog(manyDistance,manyF1,'-x')
    grid on
    xlabel('d')
    ylabel('errU1')
    title(['R_1=' num2str(rBubble1) ' R_2=' num2str(rBubble2) ' n=' num2str(manyN(1))])
    
%     subplot(2,1,2)
%     plot(manyDistance,errBubble2,'-x')
%     grid on
%     xlabel('d')
%     ylabel('errU2')
    
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








