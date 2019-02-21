%build geometry

function [x,y,minSize,maxSize] = buildGeometryPanelsParametric(tParametric,PARAM)

    elem = PARAM.n;

    %initialize
    x = cell(numel(elem),1);
    y = cell(numel(elem),1);
    minSize = zeros(numel(elem),1);
    maxSize = zeros(numel(elem),1);
    
    for i = 1:numel(elem)
       
        if PARAM.geometryPanel(i)==0    %straight line
            
            [x{i},y{i}] = lineEquation(tParametric{i},PARAM,i);
            
            xHere = x{i};   yHere = y{i};
            minSize(i) = sqrt((xHere(1)-xHere(2))^2+(yHere(1)-yHere(2))^2);
            maxSize(i) = 2*minSize(i);
            
        elseif PARAM.geometryPanel(i)==1    %arc
            
            rArc = PARAM.rArc(i);
            [x{i},y{i}] = circleEquation(tParametric{i},PARAM,i);
            
            theta = tParametric{i};
            minSize(i) = (theta(2)-theta(1))*rArc;
            maxSize(i) = 2*minSize(i);
            
        elseif PARAM.geometryPanel(i)==2 && PARAM.ellipseShape(i)==1    %ellipse
            
            %disp(['Block ' num2str(i) ' is initially an ellipse'])
            
            [xHere,yHere] = ellipseCartesian(linspace(0,pi,PARAM.n(i)+1),PARAM.D(i));
            
            rCircle = PARAM.rArc(i);
            x{i} = rCircle*xHere + PARAM.x0_Circle(i);
            y{i} = rCircle*yHere;
            
        elseif PARAM.geometryPanel(i)==2 && PARAM.ellipseShape(i)==2    %lenedre poli
            
            %disp(['Block ' num2str(i) ' is initially a sum of Legendre Polynomia'])
            
            [xHere,yHere] = spherePlus3LegendreCartesian(linspace(0,pi,PARAM.n(i)+1),PARAM.f{i},PARAM.whichF{i});
            
            rCircle = PARAM.rArc(i);
            x{i} = rCircle*xHere + PARAM.x0_Circle(i);
            y{i} = rCircle*yHere;
            
        else
            
            error('Not implemented')
            
        end
        
    end
  
end