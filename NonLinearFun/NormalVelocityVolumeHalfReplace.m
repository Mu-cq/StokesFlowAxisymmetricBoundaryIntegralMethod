%compute velocity normal to the interface

function [u,Vmiddle] = NormalVelocityVolumeHalfReplace(theta,r,q,lambda,capillary,V0,Replace,thetaReplace)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %the 2 ghost points
    theta1 = thetaReplace;
    theta2 = pi-thetaReplace;
    
    %flip
    a = [a(1:end-1) -flip(a)];
    b = [b(1:end-1) flip(b)];
    
    %compute rho in whatever position
    fVolume = @(rho) ModifyVolume3(a,b,rho,V0,Replace,thetaReplace);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    %nnn = q+1;
    xxx = [a(1:Replace-1) rMiddle.*cos(theta1) a(Replace:end-Replace+1) rMiddle.*cos(theta2) a(end-Replace+2:end)];
    yyy = [b(1:Replace-1) rMiddle.*sin(theta1) b(Replace:end-Replace+1) rMiddle.*sin(theta2) b(end-Replace+2:end)];
    
    %compute solution
    [y,N] = bem_newton_extens(xxx,yyy,q+1,lambda,capillary);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    uuu = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %only velocity without middle point
    u = [uuu(1:Replace-1); uuu(Replace+1:end-Replace); uuu(end-Replace+2:end)];
    
    %only half
    u = u(1:floor(q/2)+1);
    Vmiddle = uuu(Replace);

end