%build ellipse radius in polar coordinates

function [xp,yp,xpp,ypp] = ellipseCartesianDerivatives(theta,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse first derivative
    xp = -a*sin(theta);
    yp = b*cos(theta);
    
    %second derivative
    xpp = -a*cos(theta);
    ypp = -b*sin(theta);

end