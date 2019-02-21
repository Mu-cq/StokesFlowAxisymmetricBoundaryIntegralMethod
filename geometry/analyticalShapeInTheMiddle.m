%compute an intermediate ellipse between two given shapes

function [xMode,yMode] = analyticalShapeInTheMiddle(xMode1,yMode1,xMode2,yMode2,PARAM)

%compute gris points
[x1,y1] = fromModesToGrid(xMode1,yMode1,PARAM);
[x2,y2] = fromModesToGrid(xMode2,yMode2,PARAM);

if PARAM.shapeEllipse==1
    
    %compute D for ellipse or for legendre
    L1 = max(x1)-min(x1);   B1 = 2*max(y1);
    D1 = (L1-B1)/(L1+B1);
    L2 = max(x2)-min(x2);   B2 = 2*max(y2);
    D2 = (L2-B2)/(L2+B2);

    %compute new ellipse
    D3 = (D1+D2)/2;
    [x,y] = ellipseCartesian(pi*PARAM.t,D3);

elseif PARAM.shapeEllipse==0
    
    %compute sphere plus 
    
    
end

[xMode,yMode] = fromGridToModes(x,y,PARAM);