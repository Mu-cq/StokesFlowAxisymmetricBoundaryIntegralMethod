%build ellipse radius in polar coordinates

function [x,y,z,a,b,theta,phi,Area] = ellipseCartesian3D(theta,phi,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %create grid
    [theta,phi] = meshgrid(theta,phi);
    theta = theta(:);   phi = phi(:);
    
    %ellipse coordinates
    x = b*cos(theta).*cos(phi);
    y = b*cos(theta).*sin(phi);
    z = a*sin(theta);
    
    %surface area
    if D>=0
        e = sqrt(1-b^2/a^2);
        Area = 2*pi*b^2*(1+a/b/e*asin(e));
    else
        e = sqrt(1-a^2/b^2);
        Area = 2*pi*b^2*(1+(1-e^2)/e*atanh(e));
    end


end