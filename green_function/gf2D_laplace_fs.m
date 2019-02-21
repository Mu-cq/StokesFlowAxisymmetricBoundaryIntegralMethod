%green function of Laplace equation in 2D unbounded case, velocities and
%velocity potential for a source placed in (x0,y0)

function [phi,u,v,nxx,nxy,nyx,nyy] = gf2D_laplace_fs(x,y,x0,y0,mass)

    %coefficient for amount of mass
    A = mass/2/pi;

    %widely used variables
    dx = x-x0;
    dy = y-y0;
    r2 = dx.*dx+dy.*dy;
    r4 = r2.*r2;
    r = sqrt(r2);

    %potential
    phi = -A*log(r);
    
    %velocities
    u = A*dx./r2;
    v = A*dy./r2;
    
    %stresses
    nxx = 2*A./r2 - 4*A*dx.*dx./r4;
    nxy = - 4*A*dx.*dy./r4;
    nyx = nxy;
    nyy = 2*A./r2 - 4*A*dy.*dy./r4;

end