%build ellipse radius in polar coordinates

function [r,a,b,theta,phi] = ellipseSpherical3D(theta,phi,D)

    %minor and major axis
    b = ((1.0-D)/(1.0+D))^(1/3);
    a = 1.0/b^2;
    
    %create grid
    [theta,phi] = meshgrid(theta,phi);
    theta = theta(:);   phi = phi(:);
    
    %ellipse coordinates
    r = (a*b)./sqrt(b^2*cos(theta).^2+a^2*sin(theta).^2);
    
%     %surface area
%     if D>=0
%         e = sqrt(1-b^2/a^2);
%         [~,E] = ellipke(e^2);
%         Area = ;
%     else
%         e = sqrt(1-a^2/b^2);
%         [~,E] = ellipke(e^2);
%         Area = ;
%     end

end