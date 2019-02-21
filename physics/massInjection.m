%extensional flow

function [u,v] = massInjection(x,y,Q)

    Umass = Q./(4*pi*(x.^2+y.^2));
    theta = atan(y./x);
    theta = theta + pi*(theta<0);
    u = Umass.*cos(theta);
    v = Umass.*sin(theta);
 
end