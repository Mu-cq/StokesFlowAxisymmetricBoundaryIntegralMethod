%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralXYdns(x,y,nx,ny,DN,V0,SPECTRAL)
    
    %displace in the normal direction
    x = x + DN*nx;
    y = y + DN*ny;
    
    %dealiasing
    [x,y] = dealiasingGridXY(x,y,SPECTRAL);
    
    %residual
    R = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end