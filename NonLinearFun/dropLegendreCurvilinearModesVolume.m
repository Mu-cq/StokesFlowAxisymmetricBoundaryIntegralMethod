%drop in extensional flow: find modes of velocity in curvilnear coordinates

function [uvMode,res] = dropLegendreCurvilinearModesVolume(t,xyMode,V0,PARAM)

%legendre polyninomis
PPP = PARAM.PPP;

%allocate memory
uvMode = zeros(numel(xyMode),1);

%get shape from input
xMode = [xyMode(1); xyMode(2:2:end)];
yMode = xyMode(3:2:end);

%compute modes first y mode from volume (x is set to zero)
fVolume = @(rho) ModifyVolumeSpectralXYmodes(xMode,yMode,rho,V0,PARAM);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
yFirst = fsolve(fVolume,1,options);

yMode = [yFirst; yMode];

%get coordinates from modes
x = LegendreBuildXY(xMode,PPP);
y = LegendreBuildXY(yMode,PPP);

%compute solution
[sol,nx,ny] = bemLegendreXY(x,y,PARAM);

if PARAM.Unormal==1

    %NORMAL interaface velocities
    uNormal = sol(1:2:end-1).*nx + sol(2:2:end).*ny;
    ux = uNormal.*nx;  uy = uNormal.*ny;

elseif PARAM.Unormal==0
    
    %interaface velocities
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    
end

%compute residuals
res = sol(1:2:end-1).*nx + sol(2:2:end).*ny;
res = max(abs(res));
    
%get modes of velocities
uMode = LegendreSerieSpectralXY(ux,PPP,PARAM);
vMode = LegendreSerieSpectralXY(uy,PPP,PARAM);
    
%de-aliasing
uMode(PARAM.dealiasing+1:end) = 0;
vMode(PARAM.dealiasing+1:end) = 0;
   
%output
uvMode([1 2:2:end-1]) = uMode;
uvMode(3:2:end) = vMode(2:end);











