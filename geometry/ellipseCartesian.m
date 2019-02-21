%build ellipse radius in polar coordinates

function [x,y,L,a,b,Area] = ellipseCartesian(theta,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %ellipse coordinates
    x = a*cos(theta);
    y = b*sin(theta);
    
    %curve lenght and surface area
    if D>=0
        e = sqrt(1-b^2/a^2);
        [~,E] = ellipke(e^2);
        L = 2*a*E;
        Area = 2*pi*b^2*(1+a/b/e*asin(e));
    else
        e = sqrt(1-a^2/b^2);
        [~,E] = ellipke(e^2);
        L = 2*b*E;
        Area = 2*pi*b^2*(1+(1-e^2)/e*atanh(e));
    end

end