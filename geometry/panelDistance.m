%compute distance between two panels

function [distance,indNode,distanceNode,directionVector] = panelDistance(x,y,panel1,panel2,PARAM)

directionVector = [];

%find closest sphere point to wall
dist = distDropWall(x{panel1},y{panel1},x{panel2},y{panel2});
[distanceNode,indNode] = min(dist);

if PARAM.geometryPanel(panel1)==0 % is a straight line
    
    xLine = x{panel1};
    yLine = y{panel1};
    
    if PARAM.geometryPanel(panel2)==0 % is a straight line
        
        %in this case compute it using distance bwtween nodes
        dist = distWallDrop(x{panel1},y{panel1},x{panel2},y{panel2});
        distance = min(dist);
        
    elseif PARAM.geometryPanel(panel2)==1 % is an arc
        
        xcmBubble = PARAM.x0_Circle(panel2);
        ycmBubble = PARAM.y0_Circle(panel2);
        rBubble = PARAM.rArc(panel2);
        
        %compute minimum distance between line and sphere
        mLine = (yLine(end)-yLine(1))/(xLine(end)-xLine(1));
        aLine = mLine;  bLine = -1; cLine = yLine(1)-xLine(1)*mLine;
        distance = abs(aLine*xcmBubble+bLine*ycmBubble+cLine)/sqrt(aLine^2+bLine^2)-rBubble;
        
        %compute intersection between line and line normal line passing
        %trough circle center
        mNormal = -1/mLine;
        xI = (ycmBubble-yLine(1)+mLine*xLine(1)-mNormal*xcmBubble)/(mLine-mNormal);
        yI = yLine(1)+mLine*(xI-xLine(1));
        
        %compute vector along the distance
        dx = xI-xcmBubble;
        dy = yI-ycmBubble;
        nx = dx/sqrt(dx^2+dy^2);
        ny = dy/sqrt(dx^2+dy^2);
        directionVector = -[nx ny];
        
    else
        
        error('Not implemented')
        
    end
    
elseif PARAM.geometryPanel(panel1)==1 % is an arc
    
    xcmBubble = PARAM.x0_Circle(panel1);
    ycmBubble = PARAM.y0_Circle(panel1);
    rBubble = PARAM.rArc(panel1);
    
    if PARAM.geometryPanel(panel2)==0 % is a straight line
        
        xLine = x{panel2};   yLine = y{panel2};
        
        %compute minimum distance between line and sphere
        mLine = (yLine(end)-yLine(1))/(xLine(end)-xLine(1));
        aLine = mLine;  bLine = -1; cLine = yLine(1)-xLine(1)*mLine;
        distance = abs(aLine*xcmBubble+bLine*ycmBubble+cLine)/sqrt(aLine^2+bLine^2)-rBubble;
        
        %compute intersection between line and line normal line passing
        %trough circle center
        mNormal = -1/mLine;
        xI = (ycmBubble-yLine(1)+mLine*xLine(1)-mNormal*xcmBubble)/(mLine-mNormal);
        yI = yLine(1)+mLine*(xI-xLine(1));
        
        %compute vector along the distance
        dx = xcmBubble-xI;
        dy = ycmBubble-yI;
        nx = dx/sqrt(dx^2+dy^2);
        ny = dy/sqrt(dx^2+dy^2);
        directionVector = -[nx ny];
        
    elseif PARAM.geometryPanel(panel2)==1 % is an arc
        
        xcmBubble2 = PARAM.x0_Circle(panel2);
        ycmBubble2 = PARAM.y0_Circle(panel2);
        rBubble2 = PARAM.rArc(panel2);
        
        %compute minimum distance between sphere and sphere
        distance = sqrt((xcmBubble2-xcmBubble)^2+(ycmBubble2-ycmBubble)^2)-rBubble-rBubble2;       
        
        %compute vector along the distance
        dx = xcmBubble-xcmBubble2;
        dy = ycmBubble-ycmBubble2;
        m = dy/dx;
        thetaNormal = atan(m);
        if dx<0 && dy>0
            thetaNormal = thetaNormal+pi;
        elseif dx<0 && dy==0
            thetaNormal = pi;
        end
        directionVector = -[cos(thetaNormal) sin(thetaNormal)];
        
    else
        
        error('Not implemented')
        
    end
    
else
    
    error('Not implemented')
    
end