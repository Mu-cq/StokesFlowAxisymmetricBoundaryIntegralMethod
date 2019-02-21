%compute velocity normal to the interface

function uDropFrame = NormalVelocitySpectralVolumeFrame2(theta,r,V0,CM,Replace1,thetaReplace1,Replace2,thetaReplace2,PARAM)

    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolumeSpectralCM(theta,r,rho,V0,CM,Replace1,thetaReplace1,Replace2,thetaReplace2,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    rTwo = fsolve(fVolume,[1 1],options);
    r1 = rTwo(1);   r2 = rTwo(2);
    
    r = [r(1:Replace1-1); r1; r(Replace1:Replace2-2); r2; r(Replace2-1:end)];
    theta = [theta(1:Replace1-1); thetaReplace1; theta(Replace1:Replace2-2); thetaReplace2; theta(Replace2-1:end)];
    
    %compute solution
    [y,N] = bemSpectralTheta(theta,r,PARAM);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(:,1).*ux + N(:,2).*uy;

    %compute velocity of the droplet
    Vdrop = DropVelocityAxisSpectralTheta(r,theta,u,PARAM);
    
    %norrmal velocity in droplet frame
    uDropFrame = N(:,1).*(ux-Vdrop) + N(:,2).*uy;
    
    %take only velocity that I care about
    uDropFrame = [uDropFrame(1:Replace1-1); uDropFrame(Replace1+1:Replace2-1); uDropFrame(Replace2+1:end)];

end