%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolume3(a,b,rho,V0,Replace,theta)

    theta1 = theta;
    theta2 = pi-theta;

    x = [a(1:Replace-1) rho.*cos(theta1) a(Replace:end-Replace-1) rho.*cos(theta2) a(end-Replace:end)];
    y = [b(1:Replace-1) rho.*sin(theta1) b(Replace:end-Replace-1) rho.*sin(theta2) b(end-Replace:end)];
    
    R = axis_int_gauss_vect(x,y)-V0;

end