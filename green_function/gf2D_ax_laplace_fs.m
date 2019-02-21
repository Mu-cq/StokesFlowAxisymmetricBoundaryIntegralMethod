%green function of Laplace equation in 3D unbounded case, velocities and
%velocity potential for a source placed in (x0,0)

function [u,v,nxx,nxy,nyx,nyy] = gf2D_ax_laplace_fs(x,y,x0,mass)

    %coefficient for amount of mass
    A = mass/4/pi;

    %widely used variables
    dx = x-x0;
    dy = y;
    dx2 = dx.*dx;
    dy2 = dy.*dy;
    r = sqrt(dx.*dx+dy.*dy);
    r3 = r.^3;
    r5 = r.^5;
    
    %velocities
    u = A*dx./r3;
    v = A*dy./r3;
    
    %stresses
    nxx = 2*A*(1./r3-3*dx2./r5);
    nxy = -6*A*dx.*dy./r5;
    nyx = nxy;
    nyy = 2*A*(1./r3-3*dy2./r5);

end