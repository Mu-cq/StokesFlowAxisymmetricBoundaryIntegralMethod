%split elements when one elemnt is close to the other

function tParametric = remeshSplitProximity3(x,y,tParametric,PARAM,block1,block2)

%add nodes
if PARAM.geometryPanel(block1)==0        %remesh line
    
    %error('Not implemented')
    
    fLine = @(unk) lineEquation(unk,PARAM,block1);
    tHere = tParametric{block1};
    [xHere,yHere] = fLine(tHere);
    tMiddle = (tHere(1:end-1)+tHere(2:end))/2;
    [xMiddle,yMiddle] = fLine(tMiddle);
    
    %compute size of the element
    sizeElem = sqrt(diff(xHere).^2+diff(yHere).^2);
    
    if PARAM.geometryPanel(block2)==0    %other panel is a line
        
        error('not implemented')
        
        %compute distances between middle point and the cricle
        distance = abs(a*xMiddle+b*yMiddle+c)/sqrt(a^2+b^2);
        
        %find elemenets which are closer to the treshold
        ind1 = find(distance<sizeElem*PARAM.coeffDist);
        
        %find elemenets which are larger than the maximum allowed
        ind2 = find(sizeElem>PARAM.maxElem(block1));
        
        %new parametric coordinates
        tAdd1 = tMiddle(ind1);
        tAdd2 = tMiddle(ind2);
        tParametric{block1} = sort([tHere tAdd1 tAdd2]);
        
        %display(['Add ' num2str(numel(ind1)+numel(ind2)) ' elements to the bubble'])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fCircle(tParametric{block1});
        
    elseif PARAM.geometryPanel(block2)==1    %other panel is an arc
        
        %compute analytical formula for the panel
        %[xCircle,yCircle] = circleEquation(linspace(PARAM.thetaStart(block2),PARAM.thetaEnd(block2),1000),PARAM,block2);
        x0_circle = PARAM.x0_Circle(block2);
        y0_circle = PARAM.y0_Circle(block2);
        rArcCircle = PARAM.rArc(block2);
        
        %compute distance
        %distance = distDropWall(xMiddle,yMiddle,xCircle,yCircle);
        distance = abs(sqrt((x0_circle-xMiddle).^2+(y0_circle-yMiddle).^2)-rArcCircle);
        
        %find elemenets which are closer to the treshold
        ind1 = find(distance<sizeElem*PARAM.coeffDist);
        
        %find elemenets which are larger than the maximum allowed
        ind2 = find(sizeElem>PARAM.maxElem(block1));
        
        %new parametric coordinates
        tAdd1 = tMiddle(ind1);
        tAdd2 = tMiddle(ind2);
        tParametric{block1} = sort([tHere tAdd1 tAdd2]);
        
        %display(['Add ' num2str(numel(ind1)+numel(ind2)) ' elements to the bubble'])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fLine(tParametric{block1});
        
    else
        
        error('Not implemented')
        
    end
    
elseif PARAM.geometryPanel(block1)==1    %remesh circle
    
    %compute middle point of the panels
    fCircle = @(unk) circleEquation(unk,PARAM,block1);
    tHere = tParametric{block1};
    tMiddle = (tHere(1:end-1)+tHere(2:end))/2;
    [xMiddle,yMiddle] = fCircle(tMiddle);
    
    %compute size of the element
    sizeElem = diff(tHere)*PARAM.rArc(block1);
    
    if PARAM.geometryPanel(block2)==0    %other panel is a line
        
        %compute analytical formula for the panel
        [a,b,c] = lineCoeff(PARAM,block2);
        
        %compute distances between middle point and straight line
        distance = abs(a*xMiddle+b*yMiddle+c)/sqrt(a^2+b^2);
        
        %find elemenets which are closer to the treshold
        ind1 = find(distance<sizeElem*PARAM.coeffDist);
        
        %find elemenets which are larger than the maximum allowed
        ind2 = find(sizeElem>PARAM.maxElem(block1));
        
        %new parametric coordinates
        tAdd1 = tMiddle(ind1);
        tAdd2 = tMiddle(ind2);
        tParametric{block1} = sort([tHere tAdd1 tAdd2]);
        
        %display(['Add ' num2str(numel(ind1)+numel(ind2)) ' elements to the bubble'])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fCircle(tParametric{block1});
        
    elseif PARAM.geometryPanel(block2)==1    %other panel is a circle
        
        %compute analytical formula for the panel
        x0_circle = PARAM.x0_Circle(block2);
        y0_circle = PARAM.y0_Circle(block2);
        rArcCircle = PARAM.rArc(block2);
        
        %compute distance
        distance = abs(sqrt((x0_circle-xMiddle).^2+(y0_circle-yMiddle).^2)-rArcCircle);
        
        %find elemenets which are closer to the treshold
        ind1 = find(distance<sizeElem*PARAM.coeffDist);
        
        %find elemenets which are larger than the maximum allowed
        ind2 = find(sizeElem>2*PARAM.maxElem(block1));
        
        %new parametric coordinates
        tAdd1 = tMiddle(ind1);
        tAdd2 = tMiddle(ind2);
        tParametric{block1} = sort([tHere tAdd1 tAdd2]);
        
        %display(['Add ' num2str(numel(ind1)+numel(ind2)) ' elements to the panel' num2str(block1)])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fCircle(tParametric{block1});
        
    else
        
        error('Not implemented')
        
    end
    
end

%remove nodes
if PARAM.geometryPanel(block1)==0        %remesh line
    
    fLine = @(unk) lineEquation(unk,PARAM,block1);
    tHere = tParametric{block1};
    [xHere,yHere] = fLine(tHere);
    
    %compute size of the element
    sizeElem = sqrt(diff(xHere).^2+diff(yHere).^2);
    
    if PARAM.geometryPanel(block2)==1    %other panel is an arc
        
        %compute analytical formula for the panel
        [xCircle,yCircle] = circleEquation(linspace(PARAM.thetaStart(block2),PARAM.thetaEnd(block2),1000),PARAM,block2);
        
        %compute distance
        distance = distDropWall(xHere,yHere,xCircle,yCircle);
        
        %find elemenets which are closer to the treshold
        ind1 = ([sizeElem sizeElem(end)]>PARAM.maxElem(block1)/2);
        ind2 = (distance<[sizeElem sizeElem(end)]*PARAM.coeffDist*4);
        
        count = 1;
        for l = 1:numel(tHere)
            
           if  ind1(l)==1 || ind2(l)==1
               
              tLeave(count) = tHere(l);
               
           else
               
                count = count-1;
               
           end
           count = count+1;
            
        end
        tLeave = sort(tLeave);
        tLeave([1 end]) = [0 1];
        tParametric{block1} = tLeave;
        
        %display(['Remove ' num2str(numel(tHere)-numel(tLeave)) ' elements from the bubble'])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fLine(tParametric{block1});
        
    else
        
        error('Not implemented')
        
    end
    
elseif PARAM.geometryPanel(block1)==1    %remesh circle
    
    %compute middle point of the panels
    fCircle = @(unk) circleEquation(unk,PARAM,block1);
    tHere = tParametric{block1};
    [xHere,yHere] = fCircle(tHere);
    
    %compute size of the element
    sizeElem = diff(tHere)*PARAM.rArc(block1);
    
    if PARAM.geometryPanel(block2)==0    %other panel is a line
        
        %compute analytical formula for the panel
        [a,b,c] = lineCoeff(PARAM,block2);
        
        %compute distances between middle point and straight line
        distance = abs(a*xHere+b*yHere+c)/sqrt(a^2+b^2);
        
        %find elemenets which are closer to the treshold
        ind1 = ([sizeElem sizeElem(end)]>PARAM.minSizeElemRemesh(block1)/2);
        ind2 = (distance<[sizeElem sizeElem(end)]*PARAM.coeffDist*4);
        
        count = 1;
        for l = 1:numel(tHere)
            
           if  ind1(l)==1 || ind2(l)==1
               
              tLeave(count) = tHere(l);
               
           else
               
                count = count-1;
               
           end
           count = count+1;
            
        end
        tParametric{block1} = sort(tLeave);
        
        %display(['Remove ' num2str(numel(tHere)-numel(tLeave)) ' elements from the bubble'])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fCircle(tParametric{block1});
        
    elseif PARAM.geometryPanel(block2)==1    %other panel is a circle
        
        %compute analytical formula for the panel
        x0_circle = PARAM.x0_Circle(block2);
        y0_circle = PARAM.y0_Circle(block2);
        rArcCircle = PARAM.rArc(block2);
        
        %compute distance
        distance = abs(sqrt((x0_circle-xHere).^2+(y0_circle-yHere).^2)-rArcCircle);
        
        %find elemenets which are closer to the treshold
        ind1 = ([sizeElem sizeElem(end)]>PARAM.minSizeElemRemesh(block1)/2);
        ind2 = (distance<[sizeElem sizeElem(end)]*PARAM.coeffDist*4);
        
        count = 1;
        for l = 1:numel(tHere)
            
           if  ind1(l)==1 || ind2(l)==1
               
              tLeave(count) = tHere(l);
               
           else
               
                count = count-1;
               
           end
           count = count+1;
            
        end
        tParametric{block1} = sort(tLeave);
        
        %display(['Remove ' num2str(numel(tHere)-numel(tLeave)) ' elements from the bubble'])
        
        %new cartesian coordinates
        [x{block1},y{block1}] = fCircle(tParametric{block1});
        
    else
        
        error('Not implemented')
        
    end
    
end

addElem = 1;
countLoop = 1;
while(addElem>0)

%check if the resolution varies smoothly
tHere = tParametric{block1};
%t1 = tHere(1);  t2 = tHere(end);
dt = diff(tHere);

%look if following elemnet is too large compared to previous one
ind = find(dt(2:end)>2*dt(1:end-1));
tAdd1 = (tHere(ind+1)+tHere(ind+2))/2;

%add elements
tHere = sort([tHere tAdd1]);
dt = diff(tHere);

%look if previous elemnet is too large compared to previous one
ind = find(dt(1:end-1)>2*dt(2:end));
tAdd2 = (tHere(ind)+tHere(ind+1))/2;

%add elements
tHere = sort([tHere tAdd2]);

tParametric{block1} = tHere;

addElem = numel(tAdd1)+numel(tAdd2);
        
%display(['Add ' num2str(numel(tAdd1)+numel(tAdd2)) ' elements from the panel number' num2str(block1)])
        
if countLoop>100
    display('Too many loops')
    break
end
countLoop = countLoop+1;

%new cartesian coordinates
if PARAM.geometryPanel(block1)==0
    [x{block1},y{block1}] = fLine(tParametric{block1});
elseif PARAM.geometryPanel(block1)==1
    [x{block1},y{block1}] = fCircle(tParametric{block1});
else
        error('not implemented')
end

end        
















