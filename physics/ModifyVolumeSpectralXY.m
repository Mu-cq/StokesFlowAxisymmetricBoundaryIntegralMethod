%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralXY(perturb,x,y,rho,V0,Replace,SPECTRAL)
    
    %derivatives
    D1 = SPECTRAL.D1;
    
    %normal to ghost point
    dx = x(Replace)-x(Replace-1);   dy = y(Replace)-y(Replace-1);
    norm = sqrt(dx.^2+dy.^2);
    nNew = [dy/norm -dx/norm];
    xMiddle = (x(Replace)+x(Replace-1))/2;
    yMiddle = (y(Replace)+y(Replace-1))/2;

    %compute full shape
    x = [x(1:Replace-1); xMiddle+rho*nNew(1); x(Replace:end)];
    y = [y(1:Replace-1); yMiddle+rho*nNew(2); y(Replace:end)];
    
    %compute geomtrical derivaties
    xp = D1*x;    yp = D1*y;
  
    %compute normal vector
    h = (xp.^2+yp.^2).^(0.5);
    nx = yp./h;
    ny = -xp./h;
    
    %perturb the shape
    indPert = find(perturb~=0,1,'first');
    if indPert<Replace
        x(1:end-1) = x(1:end-1) + perturb.*nx(1:end-1);
        y(1:end-1) = y(1:end-1) + perturb.*ny(1:end-1);
    else
        x(2:end) = x(2:end) + perturb.*nx(2:end);
        y(2:end) = y(2:end) + perturb.*ny(2:end);
    end
    
    %residual
    R(1) = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end