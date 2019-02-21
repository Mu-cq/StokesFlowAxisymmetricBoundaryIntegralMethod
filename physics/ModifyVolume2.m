%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolume2(a,b,rho,V0,Replace,theta)

    x = [a(1:Replace-1) rho.*cos(theta) a(Replace:end)];
    y = [b(1:Replace-1) rho.*sin(theta) b(Replace:end)];
    
    R = axis_int_gauss_vect(x,y)-V0;

end