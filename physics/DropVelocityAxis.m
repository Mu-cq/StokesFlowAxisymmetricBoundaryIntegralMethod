%compute droplet velocity in axisymmetric domain

function Vdrop = DropVelocityAxis(x,y,uNormal)

    %center of mass
    xcm = center_mass(x,y);
    
    %function to integrate
    f = uNormal.*(x'-xcm);
    
    %volume
    Volume = axis_int_gauss(x,y);
    
    %drop velocity as surface integral
    Vdrop = int_axis_spline_symmetric(x,y,f)/Volume;

end