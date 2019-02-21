%compute center of mass of a block

function xcm = centerOfMassBlockAxis(x,y,block,PARAM)

PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

%block coordinates
[xHere,yHere,nxHere,~,panelRange] = getBlockCoordinates(x,y,PARAM,block);

%integration weights
weight = integrationOnLineWeightAxis(xHere,yHere,PARAM.orderVariableStokes(panelRange(1)),PARAM.orderGeometryStokes(panelRange(1)),PARAM.SPlinesType(panelRange(1)));

%function to integrate
f1 = 2*pi*((xHere'.^2)/2.*nxHere'.*yHere');
f2 = 2*pi*(xHere'.*nxHere'.*yHere');

%integration
num = weight*f1;
den = weight*f2;

xcm = num/den;