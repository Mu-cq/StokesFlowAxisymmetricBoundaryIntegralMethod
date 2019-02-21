%get the coordinates of the points where the variables are stored

function [x0,y0,nnn] = computeSingularityLocation(xCell,yCell,PARAM)

%get number of the singularity
numSing = 0;
nnn = zeros(numel(PARAM.n),1);
for i = 1:numel(PARAM.n)
    
    if PARAM.orderVariable(i)==0         %constant element
        
        numSing = numSing + PARAM.n(i);
        nnn(i) = PARAM.n(i);
        
    elseif PARAM.orderVariable(i)==1         %linear element
        
        numSing = numSing + PARAM.n(i)+1;
        nnn(i) = PARAM.n(i)+1;
        
    end
    
end

%compute singularity location
x0 = zeros(numSing,1);
y0 = zeros(numSing,1);
startPos = 0;
for i = 1:numel(PARAM.n)
    
    if PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==0         %constant, straight element
        
        x = xCell{i};   y = yCell{i};
        xMiddle = (x(1:end-1)+x(2:end))/2;
        yMiddle = (y(1:end-1)+y(2:end))/2;
        x0(startPos+1:nnn(i)+startPos,1) = xMiddle;
        y0(startPos+1:nnn(i)+startPos,1) = yMiddle;
        
        startPos = startPos + nnn(i);
        
    elseif PARAM.orderVariable(i)==1 && PARAM.orderGeometry(i)==0     %linear, straight element
        
        x0(startPos+1:nnn(i)+startPos,1) = xCell{i};
        y0(startPos+1:nnn(i)+startPos,1) = yCell{i};
        
        startPos = startPos + nnn(i);
        
    elseif PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==1     %constant, curved element
        
        %compute spline coeff
        if PARAM.SPlinesType==1
              [ax,bx,cx,dx,ay,by,cy,dy] = spline_natural(xCell{i},yCell{i});
        elseif PARAM.SPlinesType==2
              [ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(xCell{i},yCell{i});
        end
        
        x = @(t) ax+bx*t+cx*t^2+dx*t^3;
        y = @(t) ay+by*t+cy*t^2+dy*t^3;
        
        x0 = x(0.5)';
        y0 = y(0.5)';
        
    elseif PARAM.orderVariable(i)==1 && PARAM.orderGeometry(i)==1     %linear, curved element
        
        x0(startPos+1:nnn(i)+startPos,1) = xCell{i};
        y0(startPos+1:nnn(i)+startPos,1) = yCell{i};
        
        startPos = startPos + nnn(i);
        
    end
    
end