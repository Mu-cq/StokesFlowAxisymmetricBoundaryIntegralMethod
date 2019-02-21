%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeSpectralXYmodesNormal(perturb,xMode,yMode,V0,SPECTRAL)

    %perturb the modes
    xModePert = xMode + perturb(:,1);
    yModePert = yMode + dy;
    
    if SPECTRAL.legendre==1
        
        %get nodal values from points
        x = LegendreBuildXY(xModePert,SPECTRAL.PPP);
        y = LegendreBuildXY(yModePert,SPECTRAL.PPP);
        
    elseif SPECTRAL.legendre==0
   
        %get nodal values from points
        x = chebcoeffs2chebvals(xModePert);
        y = chebcoeffs2chebvals(yModePert);
    
    end
    
    %residual volume
    R = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;
    
    %residuals in the normal direction
    

end