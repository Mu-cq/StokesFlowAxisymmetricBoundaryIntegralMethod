%drop in extensional flow: find modes of velocity in curvilnear coordinates

function [uvMode,res,yFirst] = normalVelocityCurvilinearVolumeModes(xyMode,V0,PARAM)

%allocate memory
uvMode = zeros(numel(xyMode),1);

%get shape from input
xMode = xyMode(1:2:end-1);
yMode = xyMode(2:2:end);

%add fake modes
xMode = [xMode; zeros(PARAM.n-PARAM.dealiasing+1,1)];
yMode = [yMode; zeros(PARAM.n-PARAM.dealiasing+1,1)];

%compute modes first y mode from volume (x is set to zero)
fVolume = @(rho) ModifyVolumeSpectralXYmodes2(xMode,yMode,rho,V0,PARAM);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
yFirst = fsolve(fVolume,1,options);

xMode = [0; xMode];     yMode = [yFirst; yMode];

%get coordinates from modes
x = chebcoeffs2chebvals(xMode);
y = chebcoeffs2chebvals(yMode);

% figure(10)
% hold on
% plot(x,y)
% hold off
% axis equal
% grid on

%compute solution
here = pwd;
cd(PARAM.bem)
[sol,nx,ny] = bemChebXY(x,y,PARAM);
cd(here);

%normal velocity
uNormal = sol(1:2:end-1).*nx + sol(2:2:end).*ny;

%interaface velocities
ux = chebfun(uNormal.*nx,[0 1]);  uy = chebfun(uNormal.*ny,[0 1]);

%compute residuals
res = uNormal;
res = max(abs(res));
    
%get modes of velocities
uMode = chebcoeffs(ux);
vMode = chebcoeffs(uy);

% figure(11)
% plot(x,ux(PARAM.t))
% hold on
% plot(x,uy(PARAM.t))
% hold off
% grid on
% 
% figure(12)
% loglog(abs(uMode))
% hold on
% loglog(abs(vMode))
% hold off
% axis equal
% grid on
    
%de-aliasing
uMode = uMode(1:PARAM.dealiasing);
vMode = vMode(1:PARAM.dealiasing);
   
%output
uvMode(1:2:end-1) = uMode(2:end);
uvMode(2:2:end) = vMode(2:end);




