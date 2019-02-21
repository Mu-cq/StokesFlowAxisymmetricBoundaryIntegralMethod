%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralXYmodes2(xMode,yMode,rho,V0,SPECTRAL)
    
    %compute full shape
    xMode = [0; xMode]; yMode = [rho; yMode];
    
    if SPECTRAL.legendre==1
        
        %get nodal values from points
        x = LegendreBuildXY(xMode,SPECTRAL.PPP);
        y = LegendreBuildXY(yMode,SPECTRAL.PPP);
        
    elseif SPECTRAL.legendre==0
   
        %get nodal values from points
        x = chebcoeffs2chebvals(xMode);
        y = chebcoeffs2chebvals(yMode);
    
    end
    
    %residual
    R = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;

end