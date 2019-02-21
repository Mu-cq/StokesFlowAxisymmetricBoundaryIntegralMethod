%compute velocity normal to the interface

function [u,nx,ny,x,y] = NormalVelocitySpectralVolumeXYcont(perturb,xBase,yBase,V0,Replace,PARAM,Ca)
    
    %current capillary number
    %PARAM.Ca = Ca + perturb(end);
    PARAM.Ca = Ca;

    %add fake perturbation in ghost point (it is always zero)
    perturb = [perturb(1:Replace-1); 0; perturb(Replace:end)];
    
    %compute normal vector
    [nx,ny] = normalVectorSpectral(xBase,yBase,PARAM);
    
    %perturb respect to base shape
    x = xBase + nx.*perturb(1:end-1);
    y = yBase + ny.*perturb(1:end-1);

    %compute rho in symmetry axis
    fVolume = @(unk) ModifyVolumeSpectralXY2(x,y,nx(Replace),ny(Replace),unk,V0,Replace,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    move = fsolve(fVolume,0,options);
    
    %compute full shape
    x(Replace) = x(Replace) + nx(Replace)*move;
    y(Replace) = y(Replace) + ny(Replace)*move;
    
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
    
%     figure(11)
%     plot(x,u)
%     grid on
%     hold on
    
%     %compute drop velocity
%     Vdrop = DropVelocityAxisSpectralTheta(r,theta,u,PARAM);
%     
%     %normal velocity in droplet frame
%     ux = y(1:2:end-1);  uy = y(2:2:end);
%     u = N(:,1).*(ux-Vdrop) + N(:,2).*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];

end