%compute velocity normal to the interface

function [u,Vmiddle] = NormalVelocityVolumeHalf(theta,r,q,lambda,capillary,V0)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %flip
    a = [a -flip(a)];
    b = [b flip(b)];
    
    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolume(a,b,rho,V0);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    nnn = q+1;
    xxx = [a(1:nnn/2) 0 a(nnn/2+1:end)];
    yyy = [b(1:nnn/2) rMiddle b(nnn/2+1:end)];
    
    %compute solution
    [y,N] = bem_newton_extens(xxx,yyy,q+1,lambda,capillary);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    uuu = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %only velocity without middle point
    u = [uuu(1:nnn/2); uuu((nnn)/2+2:end)];
    
    %only half
    u = u(1:floor(q/2)+1);
    Vmiddle = uuu(nnn/2+1);

end