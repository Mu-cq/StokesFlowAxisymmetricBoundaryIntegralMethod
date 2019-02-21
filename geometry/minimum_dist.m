%compute minimum distance between many point to a vector

function dist = minimum_dist(a,b,x,y)

    %compute distance between droplet nodes and wall
    A = repmat(a,numel(x),1);
    B = repmat(b,numel(x),1);
    X = repmat(x,numel(a),1)';
    Y = repmat(y,numel(a),1)';
              
    DIST = sqrt((A-X).^2 + (B-Y).^2);
    dist = min(DIST);

end