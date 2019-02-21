%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralCM(theta,r,rho,V0,CM,Replace1,theta1,Replace2,theta2,SPECTRAL)
    
    %compute full shape
    r = [r(1:Replace1-1); rho(1); r(Replace1:Replace2-2); rho(2); r(Replace2-1:end)];
    theta = [theta(1:Replace1-1); theta1; theta(Replace1:Replace2-2); theta2; theta(Replace2-1:end)];
    
    %residual
    R(1) = VolumePolarAxisSpectral(theta,r,SPECTRAL)-V0;
    R(2) = CenterMassPolarAxisSpectral(theta,r,SPECTRAL)-CM;

end