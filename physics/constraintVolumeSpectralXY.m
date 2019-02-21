%compute function which gives the modification of the volume due to the
%motion of one point

function [DisEquality,Equality] = constraintVolumeSpectralXY(xyMode,V0,SPECTRAL)
    
    %mode
    xMode = xyMode(1:2:end-1);
    yMode = xyMode(2:2:end);
    
    %from modes to grid
    [x,y] = fromModesToGrid(xMode,yMode,SPECTRAL);
    
    %residual
    DisEquality = [];
    Equality = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end