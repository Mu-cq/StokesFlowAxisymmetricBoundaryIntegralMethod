%drop in extensional flow: find modes of velocity in curvilnear coordinates

function [uvMode,res] = dropExtensChebfunCurvilinearModesVolume(t,xyMode,V0,PARAM)

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
x = chebcoeffs2chebvals(xMode);
y = chebcoeffs2chebvals(yMode);

%compute solution
[sol,nx,ny] = bemChebXY(x,y,PARAM);

if PARAM.Unormal==1

    %NORMAL interaface velocities
    uNormal = sol(1:2:end-1).*nx + sol(2:2:end).*ny;
    ux = chebfun(uNormal.*nx,[0 1]);  uy = chebfun(uNormal.*ny,[0 1]);

elseif PARAM.Unormal==0
    
    %interaface velocities
    ux = chebfun(sol(1:2:end-1),[0 1]);  uy = chebfun(sol(2:2:end),[0 1]);
    
end

%compute residuals
res = sol(1:2:end-1).*nx + sol(2:2:end).*ny;
res = max(abs(res));
    
%get modes of velocities
uMode = chebcoeffs(ux);
vMode = chebcoeffs(uy);
    
%de-aliasing
uMode(PARAM.dealiasing+1:end) = 0;
vMode(PARAM.dealiasing+1:end) = 0;
   
%output
uvMode([1 2:2:end-1]) = uMode;
uvMode(3:2:end) = vMode(2:end);











