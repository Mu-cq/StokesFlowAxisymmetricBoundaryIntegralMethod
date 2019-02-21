%compute velocity normal to the interface for rising droplet

function uDropFrame = NormalVelocityRisingVolume(theta,r,q,lambda,capillary,PARAM,V0,CM,Replace,thetaReplace)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %rHere = 1;
    
    %compute rho in whatever position
    fVolume = @(rho) ModifyVolume2(a,b,rho,V0,Replace,thetaReplace);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    %options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15);
    rMod = fsolve(fVolume,1,options);
    r1 = rMod(1);
    
    %current full shape
    xxx = [a(1:Replace-1) r1.*cos(thetaReplace) a(Replace:end)];
    yyy = [b(1:Replace-1) r1.*sin(thetaReplace) b(Replace:end)];
    
    %recenter
    xcm = center_mass_gauss(xxx,yyy);
    xxx = xxx-xcm;
    
    %compute solution
    [y,~,~,~,N] = bem(xxx,yyy,q+2,lambda,1,capillary,PARAM.way_curv,PARAM.lealpoz,PARAM.elem,PARAM);
    
    %normal velocity IN DROPLET FRAME OF REFERENCE
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %compute velocity of the droplet
    Vdrop = DropVelocityAxis(xxx,yyy,u);
    
    %norrmal velocity in droplet frame
    uDropFrame = N(1,:)'.*(ux-Vdrop) + N(2,:)'.*uy;
    
    %take only velocity that I care about
    uDropFrame = [uDropFrame(1:Replace-1); uDropFrame(Replace+1:end)];
    
end