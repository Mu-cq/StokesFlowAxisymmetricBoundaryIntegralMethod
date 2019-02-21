%build ellipse radius in polar coordinates

function rp = ellipseFirstDerivative(theta,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse first derivative
    numP = (a^2-b^2)*sin(2*theta);
    den = b^2*cos(theta).^2+a^2*sin(theta).^2;
    denP = den.^1.5;
    rp = -0.5*(a*b)*numP./denP;
    
    %second derivative
    %numPP = -2*(a^2-b^2)*cos(2*theta).*denP + 1.5*(a^2-b^2)*sin(2*theta).*sqrt(den).*numP;
    %denPP = denP.^2;
    %rpp = 0.5*a*b*numPP./denPP;

end