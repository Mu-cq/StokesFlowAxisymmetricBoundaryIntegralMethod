%build ellipse radius in polar coordinates

function [kNorth1,kEast1,kNorth2,kEast2] = ellipseCurvature(D)

    %minor and major axis
    b = ((1-D)./(1+D)).^(1/3);
    a = 1./b.^2;
    
    %north curvature meridional plane
    kNorth1 = b./a.^2;
    
    %east curvature meridional plane
    kEast1 = a./b.^2;
    
    %north curvature aimuthal plane
    kNorth2 = 1./b;
    
    %east curvature aimuthal plane
    kEast2 = kEast1;
    
end