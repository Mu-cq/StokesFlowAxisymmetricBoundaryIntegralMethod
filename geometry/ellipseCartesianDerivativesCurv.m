%build ellipse radius in polar coordinates

function [xp,yp,xpp,ypp] = ellipseCartesianDerivativesCurv(t,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse first derivative
    xp = -pi*a*sin(pi*t);
    yp = pi*b*cos(pi*t);
    
    %second derivative
    xpp = -a*cos(pi*t)*pi^2;
    ypp = -b*sin(pi*t)*pi^2;

end