%build ellipse radius in polar coordinates (from wikipedia, parametric surfaces)

function K = computeCurvature3DcosTheta(xt,yt,zt,xs,ys,zs,xtt,ytt,ztt,xss,yss,zss,xts,yts,zts,nx,ny,nz)
    
    ru = [xt'; yt'; zt'];
    rv = [xs'; ys'; zs'];
        
    ruu = [xtt'; ytt'; ztt'];
    ruv = [xts'; yts'; zts'];
    rvv = [xss'; yss'; zss'];
       
    n = [nx'; ny'; nz'];
       
    %first fundamental form
    E = dot(ru,ru);
    F = dot(ru,rv);
    G = dot(rv,rv);
       
    %second fundamental form
    L = dot(ruu,n);
    M = dot(ruv,n);
    N = dot(rvv,n);
       
    %curvature
    K = -(E.*N-2*F.*M+G.*L)'./(E.*G-F.^2)';
       
end
















