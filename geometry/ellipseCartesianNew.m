%build ellipse radius in polar coordinates

function [y,L] = ellipseCartesianNew(x,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse coordinates
    y = b*sqrt(1-x.^2./a.^2);
    
    %curve lenght
    e = sqrt(1-b^2/a^2);
    [~,E] = ellipke(e);
    L = 2*a*E;

end