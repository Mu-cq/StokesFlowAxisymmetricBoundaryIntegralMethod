%pertub shape with second legendre polynomial

function [x,y] = perturbWithSecondLegendrePoly(x,y,delta,PARAM)

display('Perturb shape with second Legendre polynomial')

%legendre second poly
P2 = legendre(2,cos(pi*PARAM.t));
P2 = P2(1,:)';

%compute normal vector
[nx,ny] = normalVectorCurvilinear(x,y,PARAM.D1);

%perturb
x = x + nx.*P2*delta;
y = y + ny.*P2*delta;

%dealiasing
[x,y] = dealiasingGridXY(x,y,PARAM);

%compute rho in symmetry axis
fVolume = @(unk) ModifyVolumeSpectralXYdns(x,y,nx,ny,unk,PARAM.V0,PARAM);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
move = fsolve(fVolume,0,options);
    
%compute full shape
x = x + nx*move;
y = y + ny*move;
    
%dealiasing
[x,y] = dealiasingGridXY(x,y,PARAM);