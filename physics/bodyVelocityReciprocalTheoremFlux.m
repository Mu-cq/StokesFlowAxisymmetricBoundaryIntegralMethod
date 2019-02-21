%compute force on panel

function [Vrec,Frec] = bodyVelocityReciprocalTheoremFlux(x,y,F,U,dfx,dfy,range,PARAM)

%change variables name
PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

%compute location of the singularity
[Xsing,Ysing,nnn] = computeSingularityLocation(x,y,PARAM);

warning('Bug in this function')

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
        [nx,ny] = computeNormalVector(x{i},y{i},PARAM.orderGeometry(i));
        
        %compute flux in x and y direction
        %vSlipHereX = vSlip(startMatrix:endMatrix).*tx';
        %vSlipHereY = vSlip(startMatrix:endMatrix).*ty';
        fluxX = PARAM.BCflux{i}.*nx';
        fluxY = PARAM.BCflux{i}.*ny';
               
        %force on panel
        toInt = fluxX.*dfx(startMatrix:endMatrix) + fluxY.*dfy(startMatrix:endMatrix);
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