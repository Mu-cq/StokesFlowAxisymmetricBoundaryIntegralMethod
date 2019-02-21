%compute velocity normal to the interface

function [uMode,nx,ny,x,y] = NormalVelocitySpectralVolumeXYmodesCont(perturb,xBase,yBase,V0,xcm,PARAM,Ca)

    PARAM.Ca = Ca;

    %add zeros because of dealiasing
    xBase(PARAM.dealiasing+1:PARAM.n+1) = 0;
    yBase(PARAM.dealiasing+1:PARAM.n+1) = 0;

    %compute current grid points
    if PARAM.legendre==1||PARAM.legendre==2
        PPP = PARAM.PPP;
        xGrid = LegendreBuildXY(xBase,PPP);
        yGrid = LegendreBuildXY(yBase,PPP);
    elseif PARAM.legendre==0
        PPP = PARAM.TTT;
        xGrid = chebcoeffs2chebvals(xBase);
        yGrid = chebcoeffs2chebvals(yBase);
    end
    
    %perturbation of the physical grid
    if PARAM.dropFrame==0
            perturb = sum(repmat(perturb,1,numel(xGrid)).*PPP(2:PARAM.dealiasing,:));
    elseif PARAM.dropFrame==1
            perturb = sum(repmat(perturb,1,numel(xGrid)).*PPP(3:PARAM.dealiasing,:));
    end

    %compute normal vector
    [nx,ny] = normalVectorSpectral(xGrid,yGrid,PARAM);
    
    %perturb respect to base shape
    xGrid = xGrid + nx.*perturb';
    yGrid = yGrid + ny.*perturb';
    
    %dealiasing
    [xGrid,yGrid] = dealiasingGridXY(xGrid,yGrid,PARAM);

    if PARAM.dropFrame==0

        %compute rho in symmetry axis
        fVolume = @(unk) ModifyVolumeSpectralXYdns(xGrid,yGrid,nx,ny,unk,V0,PARAM);
        options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
        move = fsolve(fVolume,0,options);

        %compute full shape
        x = xGrid + nx*move;
        y = yGrid + ny*move;
    
    elseif PARAM.dropFrame==1
        
        %compute rho in symmetry axis
        fVolume = @(unk) ModifyVolumeCmSpectralXY(xGrid,yGrid,nx,ny,unk,V0,xcm,PARAM);
        options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
        move = fsolve(fVolume,[0 0],options);

        %compute full shape
        x = xGrid + nx*move(1).*PPP(1,:)' + nx*move(2).*PPP(2,:)';
        y = yGrid + ny*move(1).*PPP(1,:)' + ny*move(2).*PPP(2,:)';
    
    end
    
    %dealiasing
    [x,y] = dealiasingGridXY(x,y,PARAM);
    
    %compute solution
    [sol,nxGrid,nyGrid] = bemSpectralXY(x,y,PARAM);
    
    %normal velocity
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    u = nxGrid.*ux + nyGrid.*uy;
    
    %compute the modes from the grid points
    if PARAM.legendre==1||PARAM.legendre==2

           uMode = LegendreSerieSpectralXY(u,PARAM.PPP,PARAM);

    elseif PARAM.legendre==0

           u = chebfun(u,[0 1]);
           uMode = chebcoeffs(u);

    end
    %dealiasing
    if PARAM.dropFrame==0
        uMode = uMode(2:PARAM.dealiasing);
    elseif PARAM.dropFrame==1
        uMode = uMode(3:PARAM.dealiasing);
    end

end