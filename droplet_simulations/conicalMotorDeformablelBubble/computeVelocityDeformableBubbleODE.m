% compute velocity of the motor and bubble once the geomtry is known

function [U,yStokes,x,y,PARAM] = computeVelocityDeformableBubbleODE(t,var,tParametricBase,Ca,beta,PARAM)

%get position of the motor
posMotor = var(end);

%get right parameters
PARAM.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
PARAM.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;

%bubble coordinates
nNodesBubble = (numel(var)-1)/2;
xBubble = var(1:2:2*nNodesBubble-1);
yBubble = var(2:2:2*nNodesBubble);
yBubble([1 end]) = [0 0];

%compute curvilinear bubble
dl = sqrt(diff(xBubble).^2+diff(yBubble).^2);
l = [0 cumsum(dl')];
PARAM.ppx(5) = spline(l/l(end),xBubble');
PARAM.ppy(5) = spline(l/l(end),yBubble');

%remesh
if t>0
     [xBase,yBase] = buildGeometryPanelsParametric(tParametricBase,PARAM);
     xBase = xBase(1:4); xBase{5} = xBubble';
     yBase = yBase(1:4); yBase{5} = yBubble';

     tParametric = remeshPanels(xBase,yBase,tParametricBase,tParametricBase,PARAM);
     for k = 1:numel(PARAM.n)
         PARAM.n(k) = numel(tParametric{k})-1;
     end
else
    tParametric = tParametricBase;
end

%build shape
[xBlock1,yBlock1] = buildGeometryPanelsParametric(tParametric,PARAM);
xBlock1 = xBlock1(1:4);
yBlock1 = yBlock1(1:4);
x = xBlock1; x{5} = xBubble';
y = yBlock1; y{5} = yBubble';
PARAM.n(5) = nNodesBubble-1;

% figure(100)
% hold off
% plotGeometryDrop(x,y,PARAM,1)
% drawnow

%compute volume flow rate to the bubble
PARAM.Qsource(2) = volumeFlowRate(xBubble',yBubble',Ca,beta,PARAM);

%solve Stokes equation
[yStokes,~,~,~,~,nnn] = BEM_Stokes(x,y,PARAM);

%blocks velocities
Ucone = yStokes(end);
Ububble = yStokes(2*sum(nnn(1:4))+1:end-1);

%enforce symmetry condition
Ububble([2 end]) = [0 0];

%build output
U = [Ububble; Ucone];










