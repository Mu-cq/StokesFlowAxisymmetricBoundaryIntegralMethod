%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeAreaSpectralXY(x,y,nx,ny,DN,V0,Area,SPECTRAL)
    
    if SPECTRAL.legendre==1||SPECTRAL.legendre==2
       PPP = SPECTRAL.PPP;
    elseif SPECTRAL.legendre==0
       PPP = SPECTRAL.TTT;
    end
    P1 = PPP(1,:)';
    P2 = PPP(3,:)';
    
    DN1 = DN(1);    DN2 = DN(2);

    %displace in the normal direction
    x = x + DN1*nx.*P1 + DN2*nx.*P2;
    y = y + DN1*ny.*P1 + DN2*ny.*P2;
    
    %residual
    R(1) = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)-V0;
    R(2) = surfaceCurvilinearAxisSpectral(x,y,SPECTRAL)-Area;

end