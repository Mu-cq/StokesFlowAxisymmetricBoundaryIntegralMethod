%build ellipse radius in polar coordinates

function [x,y] = ellipseCartesianCurv(t,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse coordinates
    x = a*cos(pi*t);
    y = b*sin(pi*t);

end