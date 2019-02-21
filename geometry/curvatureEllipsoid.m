%build ellipse radius in polar coordinates (fromula from Wolfram)

function K = curvatureEllipsoid(theta,phi,a,b,c)

    a2 = a^2;
    b2 = b^2;
    c2 = c^2;
    
    num = a*b*c*(3*(a2+b2) + 2*c2 + (a2+b2-2*c2)*cos(2*theta)-2*(a2-b2)*cos(2*phi).*sin(theta).^2);
    
    den = 4*(a2*b2*cos(theta).^2 + c2*(b2*cos(phi).^2+a2*sin(phi).^2).*sin(theta).^2).^1.5;

    K = num./den;
    
end