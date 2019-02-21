%build ellipse radius in polar coordinates (wikipedia)

function [nx,ny,nz,den] = computeNormalVector3D(xt,yt,zt,xs,ys,zs)
        
    ru = [xt'; yt'; zt']; 
    rv = [xs'; ys'; zs'];
       
    num = cross(ru,rv);
    den = sqrt(sum(num.^2))';
    DEN = repmat(den',3,1);
    N = num./DEN;
       
    nx = N(1,:)';
    ny = N(2,:)';
    nz = N(3,:)';

end
















