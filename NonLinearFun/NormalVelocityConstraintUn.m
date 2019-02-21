%compute velocity normal to the interface

function u = NormalVelocityConstraintUn(theta,r,q,lambda,capillary,Replace,thetaReplace)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolumeConstraintUn(a,b,rho,Replace,thetaReplace,capillary,lambda,q);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    %nnn = q+1;
    xxx = [a(1:Replace-1) rMiddle.*cos(thetaReplace) a(Replace:end)];
    yyy = [b(1:Replace-1) rMiddle.*sin(thetaReplace) b(Replace:end)];
    
%     figure
%     plot(xxx,yyy,'o-')
%     axis equal
%     grid on
%     xlabel('x')
%     ylabel('y')
    
    %compute solution
    [y,N] = bem_newton_extens(xxx,yyy,q,lambda,capillary);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];

end