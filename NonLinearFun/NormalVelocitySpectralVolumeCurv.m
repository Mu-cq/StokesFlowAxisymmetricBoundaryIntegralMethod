%compute velocity normal to the interface

function u = NormalVelocitySpectralVolumeCurv(xy,V0,Replace,PARAM)

    %compute point position such that volume is conserved
    fVolume = @(xy0) ModifyVolumeSpectralCurv(xy,xy0,V0,Replace,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    Middle = fsolve(fVolume,[-0.5 0.5],options);
    
    %reconstruct shape with all points
    x = xy(1:numel(xy)/2);  y = xy(numel(xy)/2+1:end);
    x = [x(1:Replace-1); Middle(1); x(Replace:end)];
    y = [y(1:Replace-1); Middle(2); y(Replace:end)];
    
    %compute solution
    [yy,nx,ny] = bemSpectralCurvilinear(x,y,PARAM);
    
    %normal velocity
    ux = yy(1:2:end-1);  uy = yy(2:2:end);
    u = nx.*ux + ny.*uy;
    
    %only velocity without middle point
    u = [u(1:Replace-1); u(Replace+1:end)];

end