%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralAxis(theta,r,rho,V0,Replace,thetaReplace,SPECTRAL)
    
    %compute full shape
    r = [rho(1); r(1:Replace-1); rho(2); r(Replace:end); rho(3)];
    theta = [0; theta(1:Replace-1); thetaReplace; theta(Replace:end); pi];
    
    %residual
    R(1) = VolumePolarAxisSpectral(theta,r,SPECTRAL)-V0;
    
    %compute derivative
    D1 = SPECTRAL.D1;
    rp = D1*r;
    
    %symmetry condition
    %symmetry = (rp-r.*sin(theta))./(rp+r.*cos(theta));
    
    %residuals
    R(2) = rp(1);
    R(3) = rp(end);
    %R(2) = symmetry(1);
    %R(3) = symmetry(end);

end