%compute distance between droplet and wall

function dist = computeDistance(x1,y1,x2,y2)

    X1 = repmat(x1,numel(x2),1);
    Y1 = repmat(y1,numel(y2),1);
    X2 = x2;
    Y2 = y2;
              
    DIST = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    dist = min(DIST);

end