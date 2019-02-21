%build ellipse radius in polar coordinates (wikipedia)

function Surface = computeSurface3D(xt,yt,zt,xs,ys,zs,wt,ws)

    ru = [xt'; yt'; zt']; 
    rv = [xs'; ys'; zs'];
       
    num = cross(ru,rv);
    den = sqrt(sum(num.^2))';
        
    %surface from normal vector
    Surface = (wt.*ws)'*den;

end
















