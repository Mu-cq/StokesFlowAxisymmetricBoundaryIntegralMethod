%build ellipse radius in polar coordinates

function r = ellipsePolarAandB(theta,a,b)
    
    %ellipse radius
    r = (a*b)./sqrt(b^2*cos(theta).^2+a^2*sin(theta).^2);

end