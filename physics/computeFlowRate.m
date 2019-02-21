%compute force on panel

function Q = computeFlowRate(x,y,flux,panel,PARAM)

%change variables name
PARAM.typeBC = PARAM.typeBClaplace;
PARAM.orderVariable = PARAM.orderVariableLaplace;
PARAM.orderGeometry = PARAM.orderGeometryLaplace;

%compute integration weight
if PARAM.orderGeometry(panel)==0 && PARAM.orderVariable(panel)==0
    
    Ysing = (y(1:end-1)+y(2:end))/2;
    dl = sqrt(diff(x).^2+diff(y).^2);
    weight = 2*pi*Ysing.*dl;
    
else
    
    error('Not implmented')
    
end

Q = weight*flux;