%remesh close wall by adding and removing elements

function [xWall,yWall,m1,m2,m] = remeshWallAddRemove(xWall,yWall,xDrop,yDrop,xCenter,yCenter,PARAM)

    addTotal  = 0;
    removeTotal  = 0;

    for k = 1:PARAM.subLoopMesh
    
        %select left and bottom part
        xDown = xWall(1:PARAM.m1+1);    yDown = yWall(1:PARAM.m1+1);
        xLeft = xWall(PARAM.m1+1:PARAM.m1+1+PARAM.m2);
        yLeft = yWall(PARAM.m1+1:PARAM.m1+1+PARAM.m2);

        %compute distance between drop and wall
        distDown = distWallDrop(xDrop,yDrop,(xDown(1:end-1)+xDown(2:end))/2,(yDown(1:end-1)+yDown(2:end))/2);
        distLeft = distWallDrop(xDrop,yDrop,(xLeft(1:end-1)+xLeft(2:end))/2,(yLeft(1:end-1)+yLeft(2:end))/2);

        %add and remove points
        [xDown, yDown, addD, removeD] = remeshThickWallAddRemove(xDown,yDown,distDown,PARAM.minimum,PARAM);
        [xLeft, yLeft, addL, removeL] = remeshCircleLeftAddRemove(xCenter,yCenter,xLeft,yLeft,distLeft,PARAM.minimum,PARAM);

        addTotal = addTotal + addD + addL;
        removeTotal = removeTotal + removeD + removeL;
        
        m1Old = PARAM.m1;   m2Old = PARAM.m2;
        PARAM.m1 = PARAM.m1+addD-removeD;  PARAM.m2 = PARAM.m2+addL-removeL;
        PARAM.m = PARAM.m+addD-removeD+addL-removeL;

        %update wall
        xWall = [xDown(1:end-1) xLeft xWall(m1Old+m2Old+2:end)];
        yWall = [yDown(1:end-1) yLeft yWall(m1Old+m2Old+2:end)];
    
    end
    
    display(['Add ' num2str(addTotal) ' elements to the wall'])
    display(['Remove ' num2str(removeTotal) ' elements from the wall'])
    
    %output number of elements
    m1 = PARAM.m1;  m2 = PARAM.m2;  m = PARAM.m;

end