%resize the droplet normally to the interface

function xyModeOut = resizeNormallyToInterface(xyMode,PARAM)

%often used
V0 = PARAM.V0;
xMode = xyMode(1:2:end-1);
yMode = xyMode(2:2:end);

%initialize
xyModeOut = zeros(numel(xyMode),1);

%from modes to grid
[x,y] = fromModesToGrid(xMode,yMode,PARAM);

%compute normal vector
[nx,ny] = normalVectorSpectral(x,y,PARAM);

%compute rho in symmetry axis
fVolume = @(unk) ModifyVolumeSpectralXYdns(x,y,nx,ny,unk,V0,PARAM);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
move = fsolve(fVolume,0,options);
    
%compute full shape
x = x + nx*move;
y = y + ny*move;
    
%dealiasing
[x,y] = dealiasingGridXY(x,y,PARAM);

%from grid to mode
[xMode,yMode] = fromGridToModes(x,y,PARAM);
xyModeOut(1:2:end-1) = xMode;
xyModeOut(2:2:end) = yMode;