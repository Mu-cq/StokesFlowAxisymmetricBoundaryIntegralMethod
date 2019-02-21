%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectral(theta,r,rho,V0,Replace,thetaReplace,SPECTRAL)
    
    %compute full shape
    r = [r(1:Replace-1); rho; r(Replace:end)];
    theta = [theta(1:Replace-1); thetaReplace; theta(Replace:end)];
    
    %residual
    R(1) = VolumePolarAxisSpectral(theta,r,SPECTRAL)-V0;

end