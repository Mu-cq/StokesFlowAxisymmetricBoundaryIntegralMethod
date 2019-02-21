%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralCurv(xy,xy0,V0,Replace,SPECTRAL)
    
    %compute full shape
    %x = xy(1:numel(xy)/2);  y = xy(numel(xy)/2+1:end);
    %x = [x(1:Replace-1); xy0(1); x(Replace:end)];
    %y = [y(1:Replace-1); xy0(2); y(Replace:end)];
    xy = [xy(1:Replace-1); xy0(1); xy(Replace:numel(xy)/2+1+Replace-1); xy0(2); xy(numel(xy)/2+1+Replace:end)];
    
    %residual
    x = xy(1:numel(xy)/2);  y = xy(1:numel(xy)/2);
    R(1) = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end