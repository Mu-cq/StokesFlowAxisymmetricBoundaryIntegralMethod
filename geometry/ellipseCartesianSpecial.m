%build ellipse radius in polar coordinates

function [x,y,L] = ellipseCartesianSpecial(theta,D,expon)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse coordinates
    x = a*cos(theta).^expon;
    y = b*sin(theta).^expon;
    
    %curve lenght
    e = sqrt(1-b^2/a^2);
    [~,E] = ellipke(e^2);
    L = 2*a*E;

end