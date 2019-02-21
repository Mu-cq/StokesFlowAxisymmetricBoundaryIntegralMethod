%compute slip velocity once he species concentration are known

function [vSlip,l] = computeSlipVelPhoretic(x,y,conc,PARAM,panel,D1)

%change variables name
PARAM.typeBC = PARAM.typeBClaplace;
PARAM.orderVariable = PARAM.orderVariableLaplace;
PARAM.orderGeometry = PARAM.orderGeometryLaplace;

if PARAM.orderVariable(panel)==0 && PARAM.orderGeometry(panel)==0
    
    %compute arc lenght constant elem
    dl = sqrt(diff(x).^2+diff(y).^2);
    l = zeros(numel(dl)-1,1);
    l(1) = dl(1)/2;
    for i = 2:numel(x)-1
        
        l(i) = (dl(i-1)+dl(i))/2 + l(i-1);
        
    end
    
elseif PARAM.orderVariable(panel)==1 && PARAM.orderGeometry(panel)==0
    
    %compute arc lenght constant elem
    dl = sqrt(diff(x).^2+diff(y).^2);
    l = zeros(numel(dl)+1,1);
    l(1) = 0;
    for i = 2:numel(x)
        
        l(i) = dl(i-1) + l(i-1);
        
    end
    
end

%compute derivative in the arc lenght
dcdt = D1*conc;
dldt = D1*l;
vSlip = dcdt./dldt;