%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeConstraintUn(a,b,rho,Replace,theta,capillary,lambda,q)

    x = [a(1:Replace-1) rho.*cos(theta) a(Replace:end)];
    y = [b(1:Replace-1) rho.*sin(theta) b(Replace:end)];
    
    %compute velocity
    [sol,N] = bem_newton_extens(x,y,q,lambda,capillary);
    
    %normal velocity
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %compute residuals
    R = int_axis_spline_symmetric(x,y,u);

end