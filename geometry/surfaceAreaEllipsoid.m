%build ellipse radius in polar coordinates

function S = surfaceAreaEllipsoid(a,b,c)

    %surface area
    if (c>=a) && (a>=b)
        
        cosPhi = b/c;
        sinPhi = sqrt(1-cosPhi^2);
        phiEll = acos(b/c);
        k2 = c^2*(a^2-b^2)/a^2/(c^2-b^2);
        k = sqrt(k2);
        
        F = ellipticF(phiEll,k^2);
        E = ellipticE(phiEll,k^2);
        
        S = 2*pi*b^2 + 2*pi*c*a/sin(phiEll) * (E*sinPhi^2+F*cosPhi^2);
        
    elseif (a>=c) && (c>=b)
        
        cosPhi = b/a;
        sinPhi = sqrt(1-cosPhi^2);
        phiEll = acos(b/a);
        k2 = a^2*(c^2-b^2)/c^2/(a^2-b^2);
        k = sqrt(k2);
        
        F = ellipticF(phiEll,k^2);
        E = ellipticE(phiEll,k^2);
        
        S = 2*pi*b^2 + 2*pi*c*a/sin(phiEll) * (E*sinPhi^2+F*cosPhi^2);
        
    elseif (a>=b) && (b>=c)

        phi = asin(sqrt(1-c^2/a^2));
        k = sqrt(a^2*(b^2-c^2)/b^2/(a^2-c^2));
        
        F = ellipticF(phi,k^2);
        E = ellipticE(phi,k^2);
        
        S = 2*pi*(c^2 + b*c^2/sqrt(a^2-c^2)*F+b*sqrt(a^2-c^2)*E);
        
    end

end