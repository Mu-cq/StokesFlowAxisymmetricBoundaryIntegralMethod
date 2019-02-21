%pertub shape with second legendre polynomial

function [x,y] = perturbWithThirdLegendrePoly(x,y,delta,PARAM)

display('Perturb shape with third Legendre polynomial')

%legendre second poly
P3 = legendre(3,cos(pi*PARAM.t));
P3 = P3(1,:)';

%compute normal vector
[nx,ny] = normalVectorCurvilinear(x,y,PARAM.D1);

%perturb
x = x + nx.*P3*delta;
y = y + ny.*P3*delta;

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