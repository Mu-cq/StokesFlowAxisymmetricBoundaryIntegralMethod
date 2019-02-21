%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%physical parameters
rBubble = 1;                    % radius of the bubble
manyGap = logspace(-3,0,50);
%manyGap = 1e-3;
rDot = 1;
D = 0;
elem = 400;

%options
plotShape = 1;

%grid for velocity field
oneline = 2;
if oneline==0
    
    meshSize = 100;
    subPoint = 0;   coeffSub = manyGap;
    %xxx = linspace(0,5*manyGap,meshSize);
    %yyy = linspace(0,80*manyGap,meshSize);
    xxx = linspace(0,0.02,meshSize);
    yyy = linspace(0,0.2,meshSize);
    xxx = linspace(0,4,meshSize);
    yyy = linspace(0,4,meshSize);
    [X0,Y0] = meshgrid(xxx,yyy);
    
elseif oneline==1
    
    rrr = logspace(0,4,100);
    theta = pi/4;
    X = rrr.*cos(theta);
    Y = rrr.*sin(theta);
    
end

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
PARAM.rArc = rBubble;               % theta starting point for the arc
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
PARAM.panelType = 0;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 0;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = 1;
PARAM.deflationConstant = 4*pi*rBubble^2;
PARAM.repulsiveForces = 0;

%print to screen
printToScreenStokes(PARAM)

%euler loop
F1 = zeros(numel(manyGap),1);
for i = 1:numel(manyGap)
    
    disp([num2str(i) ' of ' num2str(numel(manyGap))])
    
    %build geometry
    gap = manyGap(i);
    PARAM.x0_Circle = gap+rBubble;
    xx2 = linspace(0,1,PARAM.n(1)+1);
    temp = (exp(8*xx2)-1)/(exp(8)-1)*pi;
    temp = flip(diff(temp));
    tParametric{1} = [0 cumsum(temp)];
    [x,y,PARAM.minSize,PARAM.maxSize] = buildGeometryPanelsParametric(tParametric,PARAM);

    %function profile for BCs
    temp = tParametric{1};  temp = (temp(1:end-1)+temp(2:end))/2;
    if PARAM.typeBCstokes(1)==0
        PARAM.velBCaxial{1} = rDot;
        PARAM.velBCradial{1} = 0;
    elseif PARAM.typeBCstokes(1)==5
        PARAM.velBC{1} = rDot.*cos(temp);
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

    %solve Stokes equation
    [yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

    %compute force on bubble
    weight1 = integrationOnLineWeightAxis(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
    if PARAM.typeBCstokes(1)==5
        [nx1,ny1] = computeNormalVector(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
        fx1 = yStokes(1:2:2*PARAM.n(1)-1).*nx1';
    else
        fx1 = yStokes(1:2:2*PARAM.n(1)-1);
    end
    F1(i) = weight1*fx1;

end

%force from lubrication
if PARAM.typeBCstokes(1)==0
    Flubrication = -6*pi*rDot./manyGap;
    typeBC = 'rigid shell';
    myYlabel = 'F^{\lambda=0}/(a\dot{d})';
elseif PARAM.typeBCstokes(1)==5
    Flubrication = -1.5*pi*rDot./manyGap;
    typeBC = 'bubble';
    myYlabel = 'F^{\lambda \rightarrow \infty}/(a\dot{d})';
end

%plot gap
figure
if numel(manyGap)>20
    loglog(manyGap,F1,'-k')
else
    loglog(manyGap,F1,'.k','Markersize',30)
end
hold on
loglog(manyGap,Flubrication,'--')
xyLabelTex('\varepsilon',myYlabel)
grid on
legend('BEM','Lubrication','Location','Best')
title(typeBC,'interpreter','latex')

%plot error
% figure
% semilogx(manyGap,F1-Flubrication')
% xyLabelTex('\varepsilon','err')
% grid on
% title('Error','interpreter','latex')

%compute velocity field
if oneline==1||oneline==0

    %compute velocity field
    [row,col] = size(X0);
    uxGrid = zeros(row,col);
    uyGrid = zeros(row,col);
    pressure = zeros(row,col);
    for kkk = 1:row
        [X0(kkk,:),Y0(kkk,:),uxGrid(kkk,:),uyGrid(kkk,:),pressure(kkk,:)] = computeVelPressField(X0(kkk,:),Y0(kkk,:),x,y,yStokes,yStokes(end),0,PARAM,subPoint,coeffSub,1);
    end
    
end

if oneline==0
    
    figure
    plotFieldAxis(x,y,X0,Y0,uxGrid,uyGrid,1.5,0.2,1,90,PARAM)
    %axis equal
    %axis([min(min(X0)) max(max(X0)) -max(max(Y0)) max(max(Y0))])

elseif oneline==1
    
    figure
    loglog(sqrt(XsingPlot.^2+YsingPlot.^2),sqrt(ux.^2+uy.^2),'-x')
    xlabel('d')
    ylabel('|U|')
    grid on
    title('Velocity decay')
    
end

%plot stresses and velocities
if numel(manyGap)==1
    
   %sphere velocity
   disp(['Gap velocity is Ugap=' num2str(yStokes(end)-rDot)])
   
   %arclength and stresses
   l = tParametric{1};  lSing = (l(1:end-1)+l(2:end))/2;
   fn = yStokes(1:2:end-1);
   ut = yStokes(2:2:end);
  
   %plot stresses
   figure(5)
   plot(lSing,fn)
   hold on
   %plot(lSing,fy)
   xlabel('l')
   %ylabel('f_z,f_r')
   ylabel('f_n')
   grid on
    
   %plot velocity
   figure(6)
   plot(lSing,ut)
   hold on
   %plot(lSing,fy)
   xlabel('l')
   %ylabel('f_z,f_r')
   ylabel('u_t')
   grid on
   
end








