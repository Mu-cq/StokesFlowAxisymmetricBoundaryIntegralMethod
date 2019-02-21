%compute velocity normal to the interface in an extensional flow

function uNorm = NormalVelocityNorm(theta,r,q,lambda,capillary)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %compute solution
    [y,N]=bem_newton_extens(a,b,q,lambda,capillary);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %integral of the velocity on the surface
    uNorm = int_axis_spline_symmetric(a,b,u.^2);

end