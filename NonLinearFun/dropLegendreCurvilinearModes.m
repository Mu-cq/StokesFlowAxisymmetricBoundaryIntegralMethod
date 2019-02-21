%drop in extensional flow: find modes of velocity in curvilnear coordinates

function [uvMode,res,Vdrop] = dropLegendreCurvilinearModes(t,xyMode,PARAM)

%all Legendre Polynomial
PPP = PARAM.PPP;

%allocate memory
uvMode = zeros(numel(xyMode),1);

%get shape from input
xMode = xyMode(1:2:end-1);
yMode = xyMode(2:2:end);

%de-aliasing
xMode(PARAM.dealiasing+1) = 0;
yMode(PARAM.dealiasing+1) = 0;

%get coordinates from modes
x = LegendreBuildXY(xMode,PPP);
y = LegendreBuildXY(yMode,PPP);

%compute solution
[sol,nx,ny,res,~,~,Vdrop] = bemSpectralXY(x,y,PARAM);

if PARAM.Unormal==1 %move with velocity normal to the interface

    %NORMAL interaface velocities
    uNormal = sol(1:2:end-1).*nx + sol(2:2:end).*ny;
    ux = uNormal.*nx;  uy = uNormal.*ny;

elseif PARAM.Unormal==0 %move with velocity from the solution
    
    %interaface velocities
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    
end
    
%get modes of velocities
uMode = LegendreSerieSpectralXY(ux,PPP,PARAM);
vMode = LegendreSerieSpectralXY(uy,PPP,PARAM);
    
%de-aliasing
uMode(PARAM.dealiasing+1:end) = 0;
vMode(PARAM.dealiasing+1:end) = 0;
   
%output
uvMode(1:2:end-1) = uMode;
uvMode(2:2:end) = vMode;




