%drop in extensional flow: find modes of velocity in curvilnear coordinates

function [uvMode,res,Vdrop] = dropExtensChebfunCurvilinearModes(t,xyMode,PARAM)

%allocate memory
uvMode = zeros(numel(xyMode),1);

%get shape from input
xMode = xyMode(1:2:end-1);
yMode = xyMode(2:2:end);

%get coordinates from modes
x = chebcoeffs2chebvals(xMode);
y = chebcoeffs2chebvals(yMode);

%compute solution
[sol,nx,ny,res,~,~,Vdrop] = bemSpectralXY(x,y,PARAM);

if PARAM.Unormal==1
    
    %force symmetry condition
    %sol([2 end]) = [0 0];

    %NORMAL interaface velocities
    uNormal = sol(1:2:end-1).*nx + sol(2:2:end).*ny;
    ux = chebfun(uNormal.*nx,[0 1]);  uy = chebfun(uNormal.*ny,[0 1]);

elseif PARAM.Unormal==0
    
    %interaface velocities
    ux = chebfun(sol(1:2:end-1),[0 1]);  uy = chebfun(sol(2:2:end),[0 1]);
    
end
    
%get modes of velocities
uMode = chebcoeffs(ux);
vMode = chebcoeffs(uy);

%de-aliasing
uMode(PARAM.dealiasing+1:end) = 0;
vMode(PARAM.dealiasing+1:end) = 0;

%output
uvMode(1:2:end-1) = uMode;
uvMode(2:2:end) = vMode;






