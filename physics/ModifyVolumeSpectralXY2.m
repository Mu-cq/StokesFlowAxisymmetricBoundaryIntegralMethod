%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralXY2(x,y,nx,ny,move,V0,Replace,SPECTRAL)
    
    %new coordinates for ghost point
    x(Replace) = x(Replace) + nx*move;
    y(Replace) = y(Replace) + ny*move;
    
    %residual
    R = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end