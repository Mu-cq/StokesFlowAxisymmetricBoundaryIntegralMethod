%compute velocity normal to the interface

function [uMode,nx,ny,x,y,firstTwoMode,K1,K2] = NormalVelocitySpectralVolumeXYmodesRising(perturb,xBase,yBase,V0,xcm,PARAM)

    %add zeros because of dealiasing
    xBase(PARAM.dealiasing+1:PARAM.n+1) = 0;
    yBase(PARAM.dealiasing+1:PARAM.n+1) = 0;

    %compute current grid points
    if PARAM.legendre==1||PARAM.legendre==2
        PPP = PARAM.PPP;
        xGrid = LegendreBuildXY(xBase,PPP);
        yGrid = LegendreBuildXY(yBase,PPP);
        perturb = sum(repmat(perturb,1,numel(xGrid)).*PPP(3:PARAM.dealiasing,:));
    elseif PARAM.legendre==0
        PPP = PARAM.TTT;
        xGrid = chebcoeffs2chebvals(xBase);
        yGrid = chebcoeffs2chebvals(yBase);
        perturb = sum(repmat(perturb,1,numel(xGrid)).*PPP(3:PARAM.dealiasing,:));
    end

    %compute normal vector
    [nx,ny] = normalVectorSpectral(xGrid,yGrid,PARAM);
    
    %force symmetry condition
%     if PARAM.legendre==0||PARAM.legendre==2
%         nx([1 end]) = [1 -1];
%         ny([1 end]) = [0 0];
%     end
    
    %perturb respect to base shape
    xGrid = xGrid + nx.*perturb';
    yGrid = yGrid + ny.*perturb';
    
    %dealiasing
    [xGrid,yGrid] = dealiasingGridXY(xGrid,yGrid,PARAM);

    %compute rho in symmetry axis
    fVolume = @(unk) ModifyVolumeCmSpectralXY(xGrid,yGrid,nx,ny,unk,V0,xcm,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    move = fsolve(fVolume,[0 0],options);
    
    %compute full shape
    x = xGrid + nx*move(1).*PPP(1,:)' + nx*move(2).*PPP(2,:)';
    y = yGrid + ny*move(1).*PPP(1,:)' + ny*move(2).*PPP(2,:)';
    
    %dealiasing
    [x,y] = dealiasingGridXY(x,y,PARAM);
    
%     figure(10)
%     plot(x,y)
%     axis equal
%     grid on
%     hold on
    
    %compute solution
    here = pwd;
    cd(PARAM.bem)
    [sol,nxGrid,nyGrid,~,K1,K2] = bemSpectralXY(x,y,PARAM);
    %PARAM.cfunction = 1;
    %[sol,~,~,~,N] = bem(x',y',PARAM.n+1,PARAM.visc,1,PARAM.Ca,1,1,1,PARAM);
    %nxGrid = N(1,:)';   nyGrid = N(2,:)';
    cd(here)
    
    %normal velocity
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    u = nxGrid.*ux + nyGrid.*uy;
    
    %compute modes
    %compute the modes from the grid points
    if PARAM.legendre==1||PARAM.legendre==2

           uMode = LegendreSerieSpectralXY(u,PARAM.PPP,PARAM);

    elseif PARAM.legendre==0

           u = chebfun(u,[0 1]);
           uMode = chebcoeffs(u);

    end
    %dealiasing
    uMode = uMode(3:PARAM.dealiasing);
    firstTwoMode = uMode(1:2);

end