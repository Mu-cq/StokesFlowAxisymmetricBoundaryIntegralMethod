%compute velocity normal to the interface, the interface is reconstructed
%starting form modes coefficients

function [uMode,firstMode,u] = NormalVelocityVolumeModesXY(t,fMode,V0,PARAM)

    %symmetric modes
    %symmetric = PARAM.symmetric;

    %compute rho in symmetry axis
    fVolume = @(firstMode) ModifyVolumeModes(theta,fMode,firstMode,V0);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    firstMode = fsolve(fVolume,1,options);
    
    %build shape from mode
    f = [firstMode; fMode];
    r = chebcoeffs2chebvals(f);
    rCheb = chebfun(r,[0 pi]);
    r = rCheb(theta);
    %r = LegendreBuildShape(theta,f,symmetric);
    
    %aaa = (axis_int_gauss_vect(r.*cos(theta),r.*sin(theta))-V0)/V0
    
%     figure(9)
%     hold on
%     plot(r.*cos(theta),r.*sin(theta),'o-')
%     axis equal
%     grid on
%     drawnow
    
    %compute solution
    %[y,N] = bemSpectralTheta(theta,r,PARAM);
    [y,N] = bem_newton_extens(r.*cos(theta),r.*sin(theta),PARAM.n,PARAM.visc,PARAM.Ca);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
%     figure(10)
%     hold on
%     plot(theta,u,'o-')
%     grid on
%     hold off
%     drawnow
    
    %display('ciao')
    
    %compute velocity modes coeff
    uCheb = chebfun(u,[0 pi],'equi');
    uMode = chebcoeffs(uCheb);
    %uMode = LegendreSerie(theta',u,PARAM.dealiasing,symmetric,PARAM);
    
    %take minus 1 velocty coefficient
    %uMode(PARAM.dealiasing+1:end) = 0;
    if numel(uMode)<PARAM.dealiasing
        uMode = [uMode(2:end); zeros(PARAM.dealiasing-numel(uMode),1)]; 
    else
        uMode = uMode(2:PARAM.dealiasing);
    end
    %uMode = [uMode(2:end); uMode(1)];
    
%      figure(11)
%      
%      loglog(abs(uMode(2:2:end)),'o-')
%      grid on
%      hold on
%      drawnow

end