%compute distance between droplet and wall

function [dist,sign,nx] = geometryRepulsive(x1,y1,x2,y2)

    A1 = repmat(x1,numel(x2),1);    A1 = A1(:);
    B1 = repmat(y1,numel(y2),1);    B1 = B1(:);
    A2 = repmat(x2,numel(x1),1)';   A2 = A2(:);
    B2 = repmat(y2,numel(y1),1)';   B2 = B2(:);
              
    dist = sqrt((A1-A2).^2 + (B1-B2).^2);
    [dist,ind] = min(dist(:));
    
    nx = (A2(ind)-A1(ind))/dist;
    
    if nx>=0
        sign = 1;
    else
        sign = -1;
    end

end