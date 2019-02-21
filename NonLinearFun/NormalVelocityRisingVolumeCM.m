%compute velocity normal to the interface for rising droplet

function uDropFrame = NormalVelocityRisingVolumeCM(theta,r,q,lambda,capillary,PARAM,V0,CM,Replace1,thetaReplace1,Replace2,thetaReplace2)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %rHere = 1;
    
    %compute rho in whatever position
    fVolume = @(rho) ModifyVolumeCM(a,b,rho,V0,CM,Replace1,thetaReplace1,Replace2,thetaReplace2);
    %fVolume = @(rho) ModifyVolumeCM2(a,b,rho,V0,rHere,Replace1,thetaReplace1,Replace2,thetaReplace2);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    %options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15);
    rMod = fsolve(fVolume,[1 1],options);
    r1 = rMod(1);  r2 = rMod(2);
    %rMod = fsolve(fVolume,1,options);
    %r1 = rMod;  r2 = rHere;
    
    %current full shape
    xxx = [a(1:Replace1-1) r1.*cos(thetaReplace1) a(Replace1:Replace2-2) r2.*cos(thetaReplace2) a(Replace2-1:end)];
    yyy = [b(1:Replace1-1) r1.*sin(thetaReplace1) b(Replace1:Replace2-2) r2.*sin(thetaReplace2) b(Replace2-1:end)];
    
    %compute solution
    [y,~,~,~,N] = bem(xxx,yyy,q+3,lambda,1,capillary,PARAM.way_curv,PARAM.lealpoz,PARAM.elem,PARAM);
    
    %normal velocity IN DROPLET FRAME OF REFERENCE
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %compute velocity of the droplet
    Vdrop = DropVelocityAxis(xxx,yyy,u);
    
    %norrmal velocity in droplet frame
    uDropFrame = N(1,:)'.*(ux-Vdrop) + N(2,:)'.*uy;
    
    %take only velocity that I care about
    uDropFrame = [uDropFrame(1:Replace1-1); uDropFrame(Replace1+1:Replace2-1); uDropFrame(Replace2+1:end)];
    
end