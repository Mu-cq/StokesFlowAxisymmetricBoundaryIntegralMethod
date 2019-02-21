%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeC0(theta,r,c0,V0)

    %cratesina coordinates
    x = (r+c0).*cos(theta);
    y = (r+c0).*sin(theta);
    
    %compute residuals
    R = axis_int_gauss_vect(x,y)-V0;

end