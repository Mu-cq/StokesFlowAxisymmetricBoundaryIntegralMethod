%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralModes(theta,symmetric,fMode,firstMode,V0,SPECTRAL)
    
    %compute full shape
    f = [firstMode; fMode];
    r = LegendreBuildShape(theta,f(1:SPECTRAL.dealiasing+1),symmetric);
    
    %residual
    R(1) = VolumePolarAxisSpectral(theta,r,SPECTRAL)-V0;

end