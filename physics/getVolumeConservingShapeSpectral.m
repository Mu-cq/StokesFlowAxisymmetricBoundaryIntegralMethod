%get volume-conserving shape

function [x,y] = getVolumeConservingShapeSpectral(perturb,xBase,yBase,V0,PARAM)

if PARAM.legendre==1||PARAM.legendre==2
    PPP = PARAM.PPP;
elseif PARAM.legendre==0
    PPP = PARAM.TTT;
end

%perturbation of the physical grid
perturb = sum(repmat(perturb,1,numel(xBase)).*PPP(2:PARAM.dealiasing,:));

%compute normal vector
[nx,ny] = normalVectorSpectral(xBase,yBase,PARAM);
    
%perturb respect to base shape
x = xBase + nx.*perturb';
y = yBase + ny.*perturb';
    
%dealiasing
[x,y] = dealiasingGridXY(x,y,PARAM);

%compute rho in symmetry axis
fVolume = @(unk) ModifyVolumeSpectralXYdns(x,y,nx,ny,unk,V0,PARAM);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
move = fsolve(fVolume,0,options);

%compute full shape
x = x + nx*move;
y = y + ny*move;
    
%dealiasing
[x,y] = dealiasingGridXY(x,y,PARAM);