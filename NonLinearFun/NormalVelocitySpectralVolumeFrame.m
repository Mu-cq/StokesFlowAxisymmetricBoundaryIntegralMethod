%compute velocity normal to the interface

function u = NormalVelocitySpectralVolumeFrame(theta,r,V0,Replace,thetaReplace,rCM,PARAM)

    r = [rCM; r];
    theta = [0; theta];
    
    %Replace = Replace-1;

    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolumeSpectral(theta,r,rho,V0,Replace,thetaReplace,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rMiddle = fsolve(fVolume,1,options);
    
    %nnn = q+1;
    r = [r(1:Replace-1); rMiddle; r(Replace:end)];
    theta = [theta(1:Replace-1); thetaReplace; theta(Replace:end)];
    
    %compute solution
    [y,N] = bemSpectralTheta(theta,r,PARAM);
    
    %normal velocity in droplet frame
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(:,1).*ux + N(:,2).*uy;
    
    %compute drop velocity
    Vdrop = DropVelocityAxisSpectralTheta(r,theta,u,PARAM);
    
    %normal velocity in droplet frame
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(:,1).*(ux-Vdrop) + N(:,2).*uy;
    
    %only velocity without middle point
    u = [u(2:Replace-1); u(Replace+1:end)]-Vdrop;

end