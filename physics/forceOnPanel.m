%compute force on panel

function F = forceOnPanel(x,y,dfx,range,PARAM)

%change variables name
PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

%compute location of the singularity
[~,Ysing,nnn] = computeSingularityLocation(x,y,PARAM);

F = 0;
for i = range
      
    %range
    [startMatrix,endMatrix] = getSingRange(i,nnn);
      
    if PARAM.orderGeometry(i)==0 && PARAM.orderVariable(i)==0

        %integration in order to impose force free condition
        xHere = x{i};    yHere = y{i};
        dl = sqrt(diff(xHere).^2+diff(yHere).^2);
        weight = 2*pi*Ysing(startMatrix:endMatrix)'.*dl;
               
        %force on panel
        F = F + weight*dfx(startMatrix:endMatrix);

    elseif PARAM.orderGeometry(i)==1 && PARAM.orderVariable(i)==0

        error('Not Implemented')

    elseif PARAM.orderGeometry(i)==0 && PARAM.orderVariable(i)==1

        error('Not Implemented')

    elseif PARAM.orderGeometry(i)==1 && PARAM.orderVariable(i)==1

        error('Not Implemented')

    end
          
end