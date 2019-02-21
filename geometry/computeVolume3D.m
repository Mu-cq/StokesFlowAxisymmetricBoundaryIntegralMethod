%build ellipse radius in polar coordinates (wikipedia)

function Volume = computeVolume3D(z,xt,yt,xs,ys,wt,ws)
        
    IntVol = z.*(xt.*ys-xs.*yt);
    Volume = (wt.*ws)'*IntVol;

end
















