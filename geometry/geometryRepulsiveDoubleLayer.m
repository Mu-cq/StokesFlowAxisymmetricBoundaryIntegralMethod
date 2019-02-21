%compute distance between droplet and wall

function [dist,nx,ny,B1,B2,dl1,dl2] = geometryRepulsiveDoubleLayer(a1,b1,a2,b2)

    A1 = repmat(a1,numel(a2),1);
    B1 = repmat(b1,numel(b2),1);
    A2 = repmat(a2,numel(a1),1)';
    B2 = repmat(b2,numel(b1),1)';
    
    dl1 = sqrt(diff(a1).^2+diff(b1).^2);
    dl2 = sqrt(diff(a2).^2+diff(b2).^2);
    %dl1 = [dl1(1) dl1]; dl2 = [dl2(1) dl2];
    dl1 = repmat(dl1,numel(a2),1);
    dl2 = repmat(dl2,numel(a1),1);
              
    dist = sqrt((A1-A2).^2 + (B1-B2).^2);
    nx = (A2-A1)./dist;
    ny = (B2-B1)./dist;

end