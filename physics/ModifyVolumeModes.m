%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeModes(theta,fMode,firstMode,V0)
    
    %compute full shape
    f = [firstMode; fMode];
    r = chebcoeffs2chebvals(f);
    rCheb = chebfun(r,[0 pi]);
    r = rCheb(theta);
    
    %residual
    R(1) = axis_int_gauss_vect(r.*cos(theta),r.*sin(theta))-V0;

end