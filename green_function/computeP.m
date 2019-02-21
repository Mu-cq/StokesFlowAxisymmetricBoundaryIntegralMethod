%compute axysymmetric green's function p (forcing momentum equation)

function [Px,Py] = computeP(x,y,x0,y0)

    %[X,X0] = meshgrid(x,x0);
    %[Y,Y0] = meshgrid(y,y0);

    r = sqrt((x-x0).^2+(y-y0).^2);
    Px = -4*pi*y.*(x-x0)./r.^3;
    Py = -4*pi*y.*(y-y0)./r.^3;

%      r = sqrt((x-x0).^2+(y-y0).^2);
%      Px = -4*pi.*(x-x0)./r.^3;
%      Py = -4*pi.*(y-y0)./r.^3;

end