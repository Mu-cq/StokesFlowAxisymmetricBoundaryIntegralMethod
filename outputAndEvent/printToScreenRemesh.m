%print to screen information when starting the simulation

function printToScreenRemesh(PARAM)

if numel(PARAM.n)~=sum(PARAM.panels)
    error('Total panels do not correspond to panels in block')
end

disp('REMESH')

for i = 1:numel(PARAM.panels)
    
    if i==1
        sumPanels = 0;
    else
        sumPanels = sum(PARAM.panels(1:i-1));
    end
    for k = 1:PARAM.panels(i)
       
        numPan = k+sumPanels;
        if PARAM.remeshType(numPan)==0
            
            remeshType = 'not activated';
            
        elseif PARAM.remeshType(numPan)==1
            
            remeshType = 'performed by splitting the elements into two based on proximity to other panels';
        
        elseif PARAM.remeshType(numPan)==4
            
            if PARAM.distr(numPan)==0
                distrType = 'uniform';
            elseif PARAM.distr(numPan)==1
                distrType = 'curvature based';
            elseif PARAM.distr(numPan)==2
                distrType = 'upper wall-distance';
            elseif PARAM.distr(numPan)==3
                distrType = 'lateral wall-distance';
            end
            
            remeshType = ['performed with splines using ' distrType ' distribution (usually for deformable objects)'];
        
        elseif PARAM.remeshType(numPan)==5
            
            remeshType = 'performed by splitting the elements into two based on proximity to other panels (for deformable objects)';
        
        else
            
            error('Not implemented')
        
        end
        
        disp(['Remesh on panel ' num2str(k) ' of block ' num2str(i) ' is ' remeshType])
       
    end
    
end






