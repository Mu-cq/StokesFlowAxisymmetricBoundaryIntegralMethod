%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeCM(a,b,rho,V0,CM,Replace1,theta1,Replace2,theta2)

    x = [a(1:Replace1-1) rho(1).*cos(theta1) a(Replace1:Replace2-2) rho(2).*cos(theta2) a(Replace2-1:end)];
    y = [b(1:Replace1-1) rho(1).*sin(theta1) b(Replace1:Replace2-2) rho(2).*sin(theta2) b(Replace2-1:end)];
    
    R(1) = axis_int_gauss_vect(x,y)-V0;
    %R(2) = center_mass_gauss(x,y)-CM;
    R(2) = center_mass(x,y)-CM;

end