%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolume4(a,b,xEnd,V0)

    x = [a xEnd];
    y = [b 0];
    
    R = axis_int_gauss_vect(x,y)-V0;

end