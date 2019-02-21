%build ellipse radius in polar coordinates

function r = ellipsePolar(theta,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse radius
    r = (a*b)./sqrt(b^2*cos(theta).^2+a^2*sin(theta).^2);

end