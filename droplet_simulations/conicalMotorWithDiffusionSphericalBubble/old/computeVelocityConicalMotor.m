% compute velocity of the motor and bubble once the geomtry is known

function [u,v,yLaplace,yStokes] = computeVelocityConicalMotor(t,x,y,PARAM,K,beta)

%concentration on bubble surface from Henry's law
rBubble = PARAM.rArc(5);
PARAM.concBC{5} = K*(beta+2/rBubble);

%solve Laplace equation
yLaplace = BEM_Laplace(x,y,PARAM);

%compute mass flow rate to the bubble
Qmass = computeFlowRate(x{5},y{5},yLaplace(sum(PARAM.n(1:4))+1:sum(PARAM.n)),5,PARAM);

%compute associated radial displacemenet
Un = Qmass/(3*beta^rBubble^2+4*rBubble);
PARAM.velBC{5} = Un;

%solve Stokes equation
yStokes = BEM_Stokes(x,y,PARAM);

%blocks velocities
Ucone = yStokes(end-1);
Ububble = yStokes(end);

%build output
u = {Ucone Ucone Ucone Ucone [Ububble Un]};
v = {0 0 0 0 0};










