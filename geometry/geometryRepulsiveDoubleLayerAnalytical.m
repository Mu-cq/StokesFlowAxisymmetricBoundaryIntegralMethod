%compute distance between droplet and wall

function [dist,nx,ny,dl1,Y1] = geometryRepulsiveDoubleLayerAnalytical(xWall,yWall,xDrop,yDrop,theta)

    %compute distance between line and
    x0 = xWall(1);
    y0 = yWall(1);
    m = atan(theta);
    a = m;
    b = -1;
    c = y0-m*x0;
    
    %many radius
    Y1 = repmat(yWall,numel(yDrop),1);
    
    %distance between point and line
    dist = abs(a*xDrop+b*yDrop+c)/sqrt(a^2+b^2);
    dist = repmat(dist',1,numel(yWall));
              
    nx = -sin(theta);
    ny = cos(theta);
    
    dl1 = sqrt(diff(xWall).^2+diff(yWall).^2);
    dl1 = [dl1(1) dl1];
    dl1 = repmat(dl1,numel(xDrop),1);

end