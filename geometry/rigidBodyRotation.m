%solid body rotation

function [xCell,yCell,tParametric,x0,y0,xStart,xEnd,yStart,yEnd,thetaStart,thetaEnd] = rigidBodyRotation(xCell,yCell,tParametric,xC,yC,theta,x0,y0,xStart,xEnd,yStart,yEnd,thetaStart,thetaEnd,geometryType,range)

for i = range
    
    %coordinates
    x = xCell{i}-xC;   y = yCell{i}-yC;
    
    for k = 1:numel(x)
        
        %rigid body rotation
        [x(k),y(k)] = rotatePoint(x(k),y(k),theta);
    
    end
    
    %new coordinates
    xCell{i} = x + xC;
    yCell{i} = y + yC;
    
    %new coordinates center of mass of the circles
    x0(i) = x0(i) - xC;
    xStart(i) = xStart(i) - xC;
    xEnd(i) = xEnd(i) - xC;
    y0(i) = y0(i) - yC;
    yStart(i) = yStart(i) - yC;
    yEnd(i) = yEnd(i) - yC;
    [x0(i),y0(i)] = rotatePoint(x0(i),y0(i),theta);
    [xStart(i),yStart(i)] = rotatePoint(xStart(i),yStart(i),theta);
    [xEnd(i),yEnd(i)] = rotatePoint(xEnd(i),yEnd(i),theta);
    x0(i) = x0(i) + xC;
    xStart(i) = xStart(i) + xC;
    xEnd(i) = xEnd(i) + xC;
    y0(i) = y0(i) + yC;
    yStart(i) = yStart(i) + yC;
    yEnd(i) = yEnd(i) + yC;
    
    %new starting and ending theat of circles
    thetaStart = thetaStart + theta;
    thetaEnd = thetaEnd + theta;
    
    %new theta
    if geometryType(i)==1
        tParametric{i} = tParametric{i} + theta;
    end
    
end