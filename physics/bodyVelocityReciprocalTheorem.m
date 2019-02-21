%compute force on panel

function [Vrec,Frec] = bodyVelocityReciprocalTheorem(x,y,F,U,vSlip,dfx,dfy,range,PARAM)

%change variables name
PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

%compute location of the singularity
[Xsing,Ysing,nnn] = computeSingularityLocation(x,y,PARAM);

int = 0;
for i = range
      
    %range
    [startMatrix,endMatrix] = getSingRange(i,nnn);
      
    if PARAM.orderGeometry(i)==0 && PARAM.orderVariable(i)==0

        %integration in order to impose force free condition
        xHere = x{i};    yHere = y{i};
        dl = sqrt(diff(xHere).^2+diff(yHere).^2);
        weight = 2*pi*Ysing(startMatrix:endMatrix)'.*dl;
        
        %compute normal vector
        [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderVariable(i),PARAM.orderGeometry(i));
        tx = -ny;   ty = nx;
        
        %compute flux in x and y direction
        vSlipHereX = vSlip(startMatrix:endMatrix).*tx';
        vSlipHereY = vSlip(startMatrix:endMatrix).*ty';
        %vSlipHereX = PARAM.flux(i).*nx';
        %vSlipHereY = PARAM.flux(i).*ny';
               
        %force on panel
        toInt = vSlipHereX.*dfx(startMatrix:endMatrix) + vSlipHereY.*dfy(startMatrix:endMatrix);
        int = int + weight*toInt;

    elseif PARAM.orderGeometry(i)==1 && PARAM.orderVariable(i)==0

        error('Not Implemented')

    elseif PARAM.orderGeometry(i)==0 && PARAM.orderVariable(i)==1

        error('Not Implemented')

    elseif PARAM.orderGeometry(i)==1 && PARAM.orderVariable(i)==1

        error('Not Implemented')

    end
    
    %velocity from reciprocal theorem
    Vrec = -int/F;
    Frec = -int/U;
          
end