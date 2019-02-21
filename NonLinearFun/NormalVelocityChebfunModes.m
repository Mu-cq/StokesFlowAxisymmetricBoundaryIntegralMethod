%compute velocity normal to the interface, the interface is reconstructed
%starting form modes coefficients

function [uMode,firstMode,u] = NormalVelocityChebfunModes(theta,fMode,V0,PARAM)

    %symmetric modes
    symmetric = PARAM.symmetric;

    %compute rho in symmetry axis
    fVolume = @(firstMode) ModifyVolumeChebfunModes(symmetric,fMode,firstMode,V0,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    firstMode = fsolve(fVolume,1,options);
    
    %build shape from mode
    f = [firstMode; fMode];
    r = LegendreBuildShape(theta,f,symmetric);
    
    %aaa = (VolumePolarAxisSpectral(theta,r,PARAM)-V0)/V0
    
    figure(9)
    plot(r.*cos(theta),r.*sin(theta),'o-')
    axis equal
    grid on
    
    %compute solution
    [y,nx,ny] = bemMixChebLegTheta(r,PARAM);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = nx.*ux + ny.*uy;
    %u = nx.*ux + ny.*uy;
    
    %display('ciao')
    
    uCheb = chebfun(u,[0 pi]);
    
    figure(9)
    plot(theta,u,'o-')
    grid on
    
    %compute velocity modes coeff
    %uMode = LegendreSerieSpectral(theta,u,PARAM.n+1,symmetric,PARAM);
    uMode = legcoeffs(uCheb)';
    
    %take minus 1 velocty coefficient
    uMode = uMode(2:PARAM.dealiasing+1)';
    %uMode = [uMode(2:end); uMode(1)];

end