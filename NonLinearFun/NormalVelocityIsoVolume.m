%compute velocity normal to the interface

function u = NormalVelocityIsoVolume(theta,r,V0,Replace,thetaReplace,PARAM)

    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolume2(r'.*cos(theta),r'.*sin(theta),rho,V0,Replace,thetaReplace);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    %nnn = q+1;
    r = [r(1:Replace-1); rMiddle; r(Replace:end)];
    theta = [theta(1:Replace-1) thetaReplace theta(Replace:end)];
    
    %build shape
%     figure
%     plot(r'.*cos(theta),r'.*sin(theta),'k')
%     axis equal
    
    %compute solution
    [y,N] = bemIsoNewton(r'.*cos(theta),r'.*sin(theta),PARAM);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];

end