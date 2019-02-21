%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeAreaXY(x,y,nx,ny,DN,V0,Area,P2)
    
    DN1 = DN(1);    DN2 = DN(2);

    %displace in the normal direction
    x = x + DN1*nx + DN2*nx.*P2;
    y = y + DN1*ny + DN2*ny.*P2;
    
    %residual
    R(1) = axis_int_gauss_vect(x,y)-V0;
    R(2) = surf_gauss_vect(x,y)-Area;

end