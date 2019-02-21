%compute function which gives the modification of the volume due to the
%motion of one point

function AreaNorm = normAreaDifferenceSpectralXY(xyMode,Awanted,SPECTRAL)
    
    %mode
    xMode = xyMode(1:2:end-1);
    yMode = xyMode(2:2:end);
    
    %from modes to grid
    [x,y] = fromModesToGrid(xMode,yMode,SPECTRAL);
    
    %current area
    A = surfaceCurvilinearAxisSpectral(x,y,SPECTRAL);
    
    %norm of area difference
    AreaNorm = norm(A-Awanted);

end