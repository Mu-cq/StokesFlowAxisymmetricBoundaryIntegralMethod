%associated pressure to a Stokeslet and associated pressure with Stresslet
%for 2D Stokes, pressure in (X0,Y0) due to a point force or velocity in X

function [px,py,pixx,pixy,piyx,piyy] ...
    = PPI_2DStokes (X,Y,X0,Y0)

    %geometrical and often used variables
    r = sqrt((X-X0).^2+(Y-Y0).^2);
    dx = X0-X;
    dy = Y0-Y;
    
    %pressure field associated to the Stokeslet
    px = 2*dx./r.^2;
    py = 2*dy./r.^2;
    
    %pressure field associated to the Stresslet
    pixx = -2./r.^2 + 4*dx.*dx./r.^4;
    pixy = 4*dx.*dy./r.^4;
    piyx = pixy;
    piyy = -2./r.^2 + 4*dy.*dy./r.^4;

end