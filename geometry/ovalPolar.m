%build ellipse radius in polar coordinates

function r = ovalPolar(theta,D)

    %legendre function
    PPP = legendre(2,cos(theta'));
    P = PPP(1,:)';   %legendre polynomia
    
    %ellipse radius
    r = 1 + D*P;

end