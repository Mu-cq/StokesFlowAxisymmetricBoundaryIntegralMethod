%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeCM2(a,b,rho,V0,rHere,Replace1,theta1,Replace2,theta2)

    x = [a(1:Replace1-1) rho.*cos(theta1) a(Replace1:Replace2-2) rHere.*cos(theta2) a(Replace2-1:end)];
    y = [b(1:Replace1-1) rho.*sin(theta1) b(Replace1:Replace2-2) rHere.*sin(theta2) b(Replace2-1:end)];
    
    R = axis_int_gauss_vect(x,y)-V0;
    %R(2) = center_mass_gauss(x,y)-CM;
    %R(2) = center_mass(x,y)-CM;

end