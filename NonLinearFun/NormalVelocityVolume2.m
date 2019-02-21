%compute velocity normal to the interface

function u = NormalVelocityVolume2(theta,r,q,lambda,capillary,V0,Replace,thetaReplace)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolume2(a,b,rho,V0,Replace,thetaReplace);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    %nnn = q+1;
    xxx = [a(1:Replace-1) rMiddle.*cos(thetaReplace) a(Replace:end)];
    yyy = [b(1:Replace-1) rMiddle.*sin(thetaReplace) b(Replace:end)];
    
    %compute solution
    [y,N] = bem_newton_extens(xxx,yyy,q,lambda,capillary);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];

end