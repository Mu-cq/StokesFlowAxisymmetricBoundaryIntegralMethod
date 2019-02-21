%compute velocity normal to the interface

function [uMode,nx,ny] = NormalVelocitySpectralVolumeXYmodesMantain(perturb,xBase,yBase,V0,PARAM)

    %add zeros because of dealiasing
    xBase(PARAM.dealiasing+1:PARAM.n+1) = 0;
    yBase(PARAM.dealiasing+1:PARAM.n+1) = 0;

    %compute current grid points
    if PARAM.legendre==1
        PPP = PARAM.PPP;
        xGrid = LegendreBuildXY(xBase,PPP);
        yGrid = LegendreBuildXY(yBase,PPP);
        perturb = sum(repmat(perturb,1,numel(xGrid)).*PPP(2:PARAM.dealiasing,:));
    elseif PARAM.legendre==0
        TTT = PARAM.TTT;
        xGrid = chebcoeffs2chebvals(xBase);
        yGrid = chebcoeffs2chebvals(yBase);
        perturb = sum(repmat(perturb,1,numel(xGrid)).*TTT(2:PARAM.dealiasing,:));
    end

    %compute normal vector
    [nx,ny] = normalVectorSpectral(xGrid,yGrid,PARAM);
    
    %force symmetry condition
    if PARAM.legendre==0
        nx([1 end]) = [1 -1];
        ny([1 end]) = [0 0];
    end
    
    %perturb respect to base shape
    xGrid = xGrid + nx.*perturb';
    yGrid = yGrid + ny.*perturb';

    %compute rho in symmetry axis
    fVolume = @(unk) ModifyVolumeSpectralXYdns(xGrid,yGrid,nx,ny,unk,V0,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    move = fsolve(fVolume,0,options);
    
    %compute full shape
    x = xGrid + nx*move;
    y = yGrid + ny*move;
    
%     figure(10)
%     plot(x,y)
%     axis equal
%     grid on
%     hold on
    
    %compute solution
    here = pwd;
    cd(PARAM.bem)
    if PARAM.legendre==1
        [sol,nxGrid,nyGrid] = bemLegendreXY(x,y,PARAM);
    elseif PARAM.legendre==0
        [sol,nxGrid,nyGrid] = bemChebXY(x,y,PARAM);
    end
    cd(here)
    
    %normal velocity
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    u = nxGrid.*ux + nyGrid.*uy;
    
    %compute modes
    %compute the modes from the grid points
    if PARAM.legendre==1

           uMode = LegendreSerieSpectralXY(u,PARAM.PPP,PARAM);

    elseif PARAM.legendre==0

           u = chebfun(u,[0 1]);
           uMode = chebcoeffs(u);

    end
    %dealiasing
    uMode = uMode(2:PARAM.dealiasing);

end