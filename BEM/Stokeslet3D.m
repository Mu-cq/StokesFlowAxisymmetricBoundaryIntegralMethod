%3d Stokelet: compute flow field in (x,y) due to a point force in (x0,y0)

function [SXX,SXY,SXZ,SYX,SYY,SYZ,SZX,SZY,SZZ] = Stokeslet3D(x,y,z,x0,y0,z0)

    %often used
    r = sqrt((x-x0).^2+(y-y0).^2+(z-z0)^2);
    x = x-x0;
    y = y-y0;
    z = z-z0;

    SXX = 1./r + x.^2./r.^3;
    SXY = x.*y./r.^3;
    SXZ = x.*z./r.^3;
    
    SYX = SXY;
    SYY = 1./r + y.^2./r.^3;
    SYZ = y.*z./r.^3;

    SZX = SXY;
    SZY = SYZ;
    SZZ = 1./r + z.^2./r.^3;
    
end