% compute velocity of the motor and bubble once the geomtry is known

function [U,yStokes,yLaplace,Qmass5,Qmass6,x,y,PARAM] = computeVelocityConicalMotorTwoBubblesODE(t,var,tParametricBase,nPerLenght,PARAM,K,beta)

%get position of the motor and bubble and bubble radius
posMotor = var(1);
posBubble1 = var(2);
posBubble2 = var(3);
rBubble1 = var(4);
rBubble2 = var(5);

%get right parameters
PARAM.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
PARAM.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;
PARAM.x0_Circle(5) = posBubble1;
PARAM.rArc(5) = rBubble1;
PARAM.x0_Circle(6) = posBubble2;
PARAM.rArc(6) = rBubble2;

%remesh
if t>0
     nElem5 = round(PARAM.rArc(5)*pi*nPerLenght);
     nElem6 = round(PARAM.rArc(6)*pi*nPerLenght);
     if nElem5<10
         nElem5 = 10;
     end
     if nElem6<10
         nElem6 = 10;
     end
     tParametricBase{5} = linspace(0,pi,nElem5);
     tParametricBase{6} = linspace(0,pi,nElem6);
     [xBase,yBase] = buildGeometryPanelsParametric(tParametricBase,PARAM);
     tParametric = remeshPanels2(xBase,yBase,tParametricBase,tParametricBase,PARAM);
     for k = 1:numel(PARAM.n)
         PARAM.n(k) = numel(tParametric{k})-1;
     end
else
    tParametricBase{5} = linspace(0,pi,round(PARAM.rArc(5)*pi*nPerLenght));
    tParametricBase{6} = linspace(0,pi,round(PARAM.rArc(6)*pi*nPerLenght));
    PARAM.n(5) = numel(tParametricBase{5})-1;
    PARAM.n(6) = numel(tParametricBase{6})-1;
    tParametric = tParametricBase;
end
%build shape
[x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

% figure(100)
% hold off
% plotGeometryDrop(x,y,PARAM,1)
% axis([-2+posMotor 15+posMotor -8.5 8.5])
% drawnow

%concentration on bubble surface from Henry's law
PARAM.concBC{5} = K*(beta+2/rBubble1);
PARAM.concBC{6} = K*(beta+2/rBubble2);

%solve Laplace equation
yLaplace = BEM_Laplace(x,y,PARAM);

%compute mass flow rate to the bubble
Qmass5 = computeFlowRate(x{5},y{5},yLaplace(sum(PARAM.n(1:4))+1:sum(PARAM.n(1:5))),5,PARAM);

%compute mass flow rate to the bubble
Qmass6 = computeFlowRate(x{6},y{6},yLaplace(sum(PARAM.n(1:5))+1:sum(PARAM.n)),6,PARAM);

Un5 = Qmass5/(3*beta^rBubble1^2+4*rBubble1);
Un6 = Qmass6/(3*beta^rBubble2^2+4*rBubble2);
    
if t==0
        
    if Qmass6<0
            
        error('Smaller bubble is shrinking')
            
    end
        
end

%compute associated radial displacemenet
PARAM.velBC{5} = Un5;
PARAM.velBC{6} = Un6;

%solve Stokes equation
yStokes = BEM_Stokes(x,y,PARAM);

%blocks velocities
Ucone = yStokes(end-2);
Ububble5 = yStokes(end-1);
Ububble6 = yStokes(end);

%build output
U = [Ucone Ububble5 Ububble6 Un5 Un6]';










