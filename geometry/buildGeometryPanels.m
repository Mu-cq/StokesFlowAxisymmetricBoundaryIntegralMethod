%build geometry

function [x,y,tParametric,minSize,maxSize] = buildGeometryPanels(PARAM)

    error('Use parametric instaed')

    elem = PARAM.n;

    %initialize
    x = cell(numel(elem),1);
    y = cell(numel(elem),1);
    tParametric = cell(numel(elem),1);
    minSize = zeros(numel(elem),1);
    maxSize = zeros(numel(elem),1);
    
    for i = 1:numel(elem)
       
        if PARAM.geometryPanel(i)==0    %straight line
            
            xStart = PARAM.xStart(i);
            xEnd = PARAM.xEnd(i);
            yStart = PARAM.yStart(i);
            yEnd = PARAM.yEnd(i);
            
            [x{i},y{i}] = buildStraightLine(xStart,yStart,xEnd,yEnd,PARAM,i);
            
            tParametric{i} = linspace(0,1,PARAM.n(i)+1);
            xHere = x{i};   yHere = y{i};
            minSize(i) = sqrt((xHere(1)-xHere(2))^2+(yHere(1)-yHere(2))^2);
            maxSize(i) = 2*minSize(i);
            
        elseif PARAM.geometryPanel(i)==1    %arc
            
            x0 = PARAM.x0_Circle(i);
            y0 = PARAM.y0_Circle(i);
            rArc = PARAM.rArc(i);
            
            theta = linspace(PARAM.thetaStart(i),PARAM.thetaEnd(i),PARAM.n(i)+1);
            [x{i},y{i}] = buildArc(rArc,x0,y0,theta,PARAM,i);
            
            tParametric{i} = theta;
            minSize(i) = (theta(2)-theta(1))*rArc;
            maxSize(i) = 2*minSize(i);
            
        elseif PARAM.geometryPanel(i)==2    %ellipse
            
            [xHere,yHere] = ellipseCartesian(linspace(0,pi,PARAM.n(i)+1),PARAM.D(i));
            
            rCircle = PARAM.rArc(i);
            x{i} = rCircle*xHere;
            y{i} = rCircle*yHere;
            
        else
            
            error('Not implemented')
            
        end
        
    end
  
end