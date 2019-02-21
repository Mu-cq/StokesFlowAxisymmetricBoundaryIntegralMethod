%compute distance between droplet and wall

function [dist,nx,ny] = distBetweenLines(a1,b1,a2,b2)

    A1 = repmat(a1,numel(a2),1);
    B1 = repmat(b1,numel(b2),1);
    A2 = repmat(a2,numel(a1),1)';
    B2 = repmat(b2,numel(b1),1)';
              
    DIST = sqrt((A1-A2).^2 + (B1-B2).^2);
    NX = (A2-A1)./DIST;
    NY = (B2-B1)./DIST;
    [dist,ind] = min(DIST');
    
    %nx = NX(:,ind);
    %ny = NY(:,ind);
    
    nx = zeros(numel(ind),1);
    ny = zeros(numel(ind),1);
    for i = 1:numel(ind)
        nx = NX(i,ind(i));
        ny = NY(i,ind(i));
    end

end