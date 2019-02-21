%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeModes2(theta,fMode,firstMode,V0)
    
    %compute full shape
    f = [firstMode; fMode];
    
    %build shape
    r = LegendreBuildShape(theta,f,0);
    
    %residual
    R(1) = axis_int_gauss_vect(r'.*cos(theta'),r'.*sin(theta'))-V0;

end