%compute distance between droplet and wall

function dist = distWallDrop(aDrop,bDrop,aWall,bWall)

    Adrop = repmat(aDrop,numel(aWall),1);
    Bdrop = repmat(bDrop,numel(bWall),1);
    Awall = repmat(aWall,numel(aDrop),1)';
    Bwall = repmat(bWall,numel(bDrop),1)';
              
    DIST = sqrt((Adrop-Awall).^2 + (Bdrop-Bwall).^2);
    [dist,ind] = min(DIST');

end