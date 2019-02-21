%compute velocity normal to the interface

function u = NormalVelocitySpectralVolumeAxis(theta,r,V0,Replace,thetaReplace,PARAM)

    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolumeSpectralAxis(theta,r,rho,V0,Replace,thetaReplace,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rrr = fsolve(fVolume,[1 1 1],options);
    r1 = rrr(1);    r2 = rrr(2);    r3 = rrr(3);
    
    %add ghost point
    r = [r1; r(1:Replace-1); r2; r(Replace:end); r3];
    theta = [0; theta(1:Replace-1); thetaReplace; theta(Replace:end); pi];
    
    %compute solution
    [y,N] = bemSpectralTheta(theta,r,PARAM);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(:,1).*ux + N(:,2).*uy;
    
    %only velocity without middle point
    u = [u(2:Replace); u(Replace+2:end-1)];

end