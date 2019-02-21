%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralXYmodes(xMode,yMode,DN,V0,SPECTRAL)
    
    if SPECTRAL.legendre==1
        
        %get nodal values from points
        x = LegendreBuildXY(xMode,SPECTRAL.PPP);
        y = LegendreBuildXY(yMode,SPECTRAL.PPP);
        
    elseif SPECTRAL.legendre==0
   
        %get nodal values from points
        x = chebcoeffs2chebvals(xMode);
        y = chebcoeffs2chebvals(yMode);
    
    end
    
    %compute normal vector
    [nx,ny] = normalVectorSpectral(x,y,SPECTRAL);
    
    %displace in the normal direction
    x = x + DN*nx;
    y = y + DN*ny;
    
    %residual
    R = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end