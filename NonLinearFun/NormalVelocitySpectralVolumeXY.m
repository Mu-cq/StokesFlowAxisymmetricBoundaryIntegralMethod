%compute velocity normal to the interface

function [u,nx,ny] = NormalVelocitySpectralVolumeXY(perturb,x,y,V0,Replace,PARAM)

    %derivatives
    D1 = PARAM.D1;
    
    %normal to ghost point
    dx = x(Replace)-x(Replace-1);   dy = y(Replace)-y(Replace-1);
    norm = sqrt(dx.^2+dy.^2);
    nNew = [dy/norm -dx/norm];
    xMiddle = (x(Replace)+x(Replace-1))/2;
    yMiddle = (y(Replace)+y(Replace-1))/2;

    %compute rho in symmetry axis
    fVolume = @(rho) ModifyVolumeSpectralXY(perturb,x,y,rho,V0,Replace,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    xy = fsolve(fVolume,1,options);
    
    %compute full shape
    x = [x(1:Replace-1); xMiddle+xy*nNew(1); x(Replace:end)];
    y = [y(1:Replace-1); yMiddle+xy*nNew(2); y(Replace:end)];
    
    %compute geomtrical derivaties
    xp = D1*x;    yp = D1*y;
  
    %compute normal vector
    h = (xp.^2+yp.^2).^(0.5);
    nx = yp./h;
    ny = -xp./h;
    
    %perturb the shape
    indPert = find(perturb~=0,1,'first');
    if isempty(indPert)==0
        if indPert<Replace
            x(1:end-1) = x(1:end-1) + perturb.*nx(1:end-1);
            y(1:end-1) = y(1:end-1) + perturb.*ny(1:end-1);
        else
            x(2:end) = x(2:end) + perturb.*nx(2:end);
            y(2:end) = y(2:end) + perturb.*ny(2:end);
        end
    end
    
    %compute solution
    here = pwd;
    cd(PARAM.bem)
    if PARAM.legendre==1
        [sol,nx,ny] = bemLegendreXY(x,y,PARAM);
    elseif PARAM.legendre==0
        [sol,nx,ny] = bemChebXY(x,y,PARAM);
    end
    cd(here)
    
    %normal velocity
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    u = nx.*ux + ny.*uy;
    
%     %compute drop velocity
%     Vdrop = DropVelocityAxisSpectralTheta(r,theta,u,PARAM);
%     
%     %normal velocity in droplet frame
%     ux = y(1:2:end-1);  uy = y(2:2:end);
%     u = N(:,1).*(ux-Vdrop) + N(:,2).*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];
    nx = [nx(1:Replace-1); nx(Replace+1:end)];
    ny = [ny(1:Replace-1); ny(Replace+1:end)];

end