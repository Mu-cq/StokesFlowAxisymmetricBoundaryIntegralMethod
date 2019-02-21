% compute velocity of the motor and bubble once the geomtry is known

function [U,yStokes,yLaplace,Qmass,x,y,PARAM] = computeVelocityConicalMotorManyBubblesODE(t,var,tParametricBase,nPerLenght,PARAM,K,beta)

%number of bubbles
nBubble = (numel(var)-1)/2;
nWall = numel(PARAM.n)-nBubble;

%get position of the motor and bubble and bubble radius
posMotor = var(1);

%for many bubbles
posBubble = zeros(nBubble,1);
rBubble = zeros(nBubble,1);
for i = 1:nBubble
    posBubble(i) = var(i+1);
    rBubble(i) = var(1+nBubble+i);
    PARAM.x0_Circle(i+nWall) = posBubble(i);
    PARAM.rArc(i+nWall) = rBubble(i);
end

%get right parameters
PARAM.xStart(1:nWall) = PARAM.xStart(1:nWall)+posMotor;
PARAM.xEnd(1:nWall) = PARAM.xEnd(1:nWall)+posMotor;
PARAM.x0_Circle(1:nWall) = PARAM.x0_Circle(1:nWall)+posMotor;

%frid for bubbles
for i = 1:nBubble
        nNodes = round(PARAM.rArc(i+nWall)*pi*nPerLenght);
        if nNodes<10
            nNodes = 10;
        end
        tParametricBase{i+nWall} = linspace(0,pi,nNodes);
        PARAM.n(nWall+i) = nNodes-1;
end

%remesh
if t>0
    [xBase,yBase] = buildGeometryPanelsParametric(tParametricBase,PARAM);
    tParametric = remeshPanels2(xBase,yBase,tParametricBase,tParametricBase,PARAM);
    for k = 1:numel(PARAM.n)
        PARAM.n(k) = numel(tParametric{k})-1;
    end
else
    tParametric = tParametricBase;
end
%build shape
[x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

%concentration on bubble surface from Henry's law
for i = 5:nBubble+nWall
    PARAM.concBC{i} = K*(beta+2/rBubble(i-nWall));
end

%solve Laplace equation
yLaplace = BEM_Laplace(x,y,PARAM);

%compute mass flow rate to the bubble and associated radial displacemenet
Qmass = zeros(nBubble,1);
Un = zeros(nBubble,1);
for i = 5:nBubble+nWall
    Qmass(i-nWall) = computeFlowRate(x{i},y{i},yLaplace(sum(PARAM.n(1:i-1))+1:sum(PARAM.n(1:i))),i,PARAM);
    Un(i-nWall) = Qmass(i-nWall)/(3*beta^rBubble(i-nWall)^2+nWall*rBubble(i-nWall));
    PARAM.velBC{i} = Un(i-nWall);
    PARAM.stressBC{i} = 0;
end

figure(100)
hold off
PARAMplot = PARAM;
PARAMplot.velBC{6} = 1;
plotGeometryDrop(x,y,PARAMplot,1)
drawnow

%solve Stokes equation
yStokes = BEM_Stokes(x,y,PARAM);

%blocks velocities
Ucone = yStokes(end-nBubble);
Ububble = zeros(nBubble,1);
for i = 1:nBubble
    Ububble(i) = yStokes(end-nBubble+i);
end

%build output
U = [Ucone; Ububble; Un];










