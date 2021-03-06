%compute velocity normal to the interface

function u = NormalVelocitySpectralVolumeCont(theta,r,V0,Replace,thetaReplace,PARAM,Ca)

    %perturb Ca;
    CaOld = PARAM.Ca;
    PARAM.Ca = Ca;

    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolumeSpectral(theta,r,rho,V0,Replace,thetaReplace,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    %nnn = q+1;
    r = [r(1:Replace-1); rMiddle; r(Replace:end)];
    theta = [theta(1:Replace-1); thetaReplace; theta(Replace:end)];
    
    %compute solution
    %PARAM.n = PARAM.n + 1;
    [y,N] = bemSpectralTheta(theta,r,PARAM);
    %PARAM.n = PARAM.n - 1;
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(:,1).*ux + N(:,2).*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];

    %reset old Ca
    PARAM.Ca = CaOld;

end