%build ellipse radius in polar coordinates

function [x,y,z,a,b,c,theta,phi] = ellipseCartesian3Dabc(theta,phi,D1,D2)

    %axis
    b = ((1.0-D1)/(1.0+D1))^(1/3)*((1.0-D2)/(1.0+D2))^(1/3);
    a = b*(1+D1)/(1-D1);
    c = b*(1+D2)/(1-D2);
    
    %create grid
    [theta,phi] = meshgrid(theta,phi);
    theta = theta(:);   phi = phi(:);
    
    %ellipse coordinates
    x = a*sin(theta).*cos(phi);
    y = b*sin(theta).*sin(phi);
    z = c*cos(theta);

end