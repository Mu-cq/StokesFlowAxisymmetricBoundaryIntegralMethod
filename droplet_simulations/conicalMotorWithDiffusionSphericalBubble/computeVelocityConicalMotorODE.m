% compute velocity of the motor and bubble once the geomtry is known

function [U,yLaplace,yStokes,x,y,PARAM,Qmass] = computeVelocityConicalMotorODE(t,var,tParametricBase,nPerLenght,PARAM,Hcc,beta)

%get position of the motor and bubble and bubble radius
posMotor = var(1);
posBubble = var(2);
rBubble = var(3);

%get right parameters
PARAM.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
PARAM.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;
PARAM.x0_Circle(5) = posBubble;
PARAM.rArc(5) = rBubble;

%number of elements on bubble
nMinBubble = 10;
tParametricBase{5} = linspace(0,pi,round(PARAM.rArc(5)*pi*nPerLenght));
if numel(tParametricBase{5})<nMinBubble
   PARAM.n(5) = nMinBubble;
   tParametricBase{5} = linspace(0,pi,nMinBubble+1);
end

%remesh
if t>0% && sum(i==1:PARAM.remeshStep:loop)

     [xBase,yBase] = buildGeometryPanelsParametric(tParametricBase,PARAM);
     tParametric = remeshPanels(xBase,yBase,tParametricBase,tParametricBase,PARAM);
     for k = 1:numel(PARAM.n)
         PARAM.n(k) = numel(tParametric{k})-1;
     end
else
    PARAM.n(5) = numel(tParametricBase{5})-1;
    tParametric = tParametricBase;
end
%build shape
[x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

% figure(100)
% hold off
% plotGeometryDrop(x,y,PARAM,1)
% drawnow

%concentration on bubble surface from Henry's law
PARAM.concBC{5} = Hcc*(beta+2/rBubble);

%solve Laplace equation
yLaplace = BEM_Laplace(x,y,PARAM);

%compute mass flow rate to the bubble
Qmass = computeFlowRate(x{5},y{5},yLaplace(sum(PARAM.n(1:4))+1:sum(PARAM.n)),5,PARAM);

%compute velocity of inflation
Un = Qmass/(3*beta^rBubble^2+4*rBubble);

%compute associated radial displacemenet
PARAM.velBC{5} = Un;

%solve Stokes equation
yStokes = BEM_Stokes(x,y,PARAM);

%blocks velocities
Ucone = yStokes(end-1);
Ububble = yStokes(end);

%build output
U = [Ucone Ububble Un]';










