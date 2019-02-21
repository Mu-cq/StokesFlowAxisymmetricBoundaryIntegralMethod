%print to screen information when starting the simulation

function printToScreenLaplace(PARAM)

if numel(PARAM.n)~=sum(PARAM.panels)
    error('Total panels do not correspond to panels in block')
end

%LAPLACE SOLVER
display('LAPLACE SOLVER')
if PARAM.STlaplace==1
    display('Singularity treatment for Laplace equation is activated')
end

for i = 1:numel(PARAM.panels)
    
    if i==1
        sumPanels = 0;
    else
        sumPanels = sum(PARAM.panels(1:i-1));
    end
    for k = 1:PARAM.panels(i)
       
        numPan = k+sumPanels;
        if PARAM.typeBClaplace(numPan)==1
            BC = 'prescribed concentration';
        elseif PARAM.typeBClaplace(numPan)==2
            BC = 'prescribed flux';
        end
        
        if PARAM.orderGeometryLaplace(numPan)==0
            orderGeo = 'straight';
        elseif PARAM.orderGeometryLaplace(numPan)==1
            orderGeo = 'curved';
        end
        
        if PARAM.orderVariableLaplace(numPan)==0
            orderVar = 'P0';
        elseif PARAM.orderVariableLaplace(numPan)==1
            orderVar = 'P1';
        end
        
        display(['Panel ' num2str(k) ' of block ' num2str(i) ' is initialized with ' num2str(PARAM.n(numPan)) ' ' orderVar ' ' orderGeo ' elements, BC is ' BC])
       
    end
    
end

end