%compute distance between droplet and wall

function dist = computeMinimumDistance(x1,y1,x2,y2)

    X1 = repmat(x1,numel(x2),1);
    Y1 = repmat(y1,numel(y2),1);
    X2 = repmat(x2,numel(x1),1)';
    Y2 = repmat(y2,numel(y1),1)';
              
    DIST = sqrt((X1-X2).^2 + (Y1-Y2).^2);
    dist = min(min(DIST));

end