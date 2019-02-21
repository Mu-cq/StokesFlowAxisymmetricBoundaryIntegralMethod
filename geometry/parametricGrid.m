%build geometry

function tParametric = parametricGrid(PARAM)

    elem = PARAM.n;

    %initialize
    tParametric = cell(numel(elem),1);
    
    for i = 1:numel(elem)
       
        if PARAM.geometryPanel(i)==0    %straight line
            
            tParametric{i} = linspace(0,1,PARAM.n(i)+1);
            
        elseif PARAM.geometryPanel(i)==1    %arc
            
            tParametric{i} = linspace(PARAM.thetaStart(i),PARAM.thetaEnd(i),PARAM.n(i)+1);
            
        elseif PARAM.geometryPanel(i)==2    %ellipse
            
            tParametric{i} = linspace(0,1,PARAM.n(i)+1);
            
        else
            
            error('Not implemented')
            
        end
        
    end
  
end