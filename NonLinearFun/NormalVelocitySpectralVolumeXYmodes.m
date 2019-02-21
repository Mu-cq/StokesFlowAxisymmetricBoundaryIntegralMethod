%compute velocity normal to the interface

function [uMode,nx,ny,x,y] = NormalVelocitySpectralVolumeXYmodes(perturb,xMode,yMode,V0,PARAM)

    if PARAM.legendre==1
        PPP = PARAM.PPP;
    end
    
    %dealiasing
    xMode(PARAM.dealiasing+1:PARAM.n+1) = 0;
    yMode(PARAM.dealiasing+1:PARAM.n+1) = 0;
    
    %compute modes of normal vector for starting shape
    [nx,ny] = normalVectorSpectralModes(xMode,yMode,PARAM);
    
    %dealiasing
    nx = nx(2:PARAM.dealiasing);
    ny = ny(2:PARAM.dealiasing);
    
    %current shape
    xMode(2:PARAM.dealiasing) = xMode(2:PARAM.dealiasing) + perturb.*nx;
    yMode(2:PARAM.dealiasing) = yMode(2:PARAM.dealiasing) + perturb.*ny;
    xMode = xMode(2:end);
    yMode = yMode(2:end);

    %adjust first y mode for volume consevation
    fVolume = @(unk) ModifyVolumeSpectralXYmodes([0; xMode],yMode,unk,V0,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','iter');
    yFirst = fsolve(fVolume,1,options);
    
    %all modes
    xMode = [0; xMode];
    yMode = [yFirst; yMode];
    
    %full shape
    if PARAM.legendre==0
        x = chebcoeffs2chebvals(xMode);
        y = chebcoeffs2chebvals(yMode);
    elseif PARAM.legendre==1
        x = LegendreBuildXY(xMode,PPP);
        y = LegendreBuildXY(yMode,PPP);
    end
    
    figure(10)
    plot(x,y)
    axis equal
    grid on
    hold on
    
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
    
    %compute modes of normal velocity
    if PARAM.legendre==0
        u = chebfun(u,[0 1]);
        uMode = chebcoeffs(u);
        u = u(PARAM.t);
    elseif PARAM.legendre==1
        uMode = LegendreSerieSpectralXY(u,PPP,PARAM);
    end
    uMode = uMode(2:PARAM.dealiasing);
    
    figure(11)
    plot(x,u)
    grid on
    hold on
    
%     %compute drop velocity
%     Vdrop = DropVelocityAxisSpectralTheta(r,theta,u,PARAM);
%     
%     %normal velocity in droplet frame
%     ux = y(1:2:end-1);  uy = y(2:2:end);
%     u = N(:,1).*(ux-Vdrop) + N(:,2).*uy;

end