%pertub shape with second legendre polynomial

function [x,y] = perturbWithLegendrePoly(x,y,delta,PARAM)

display(['Perturb shape with ' num2str(PARAM.whichmode) 'th Legendre polynomial, amplitude=' num2str(delta)])

%legendre second poly
P = legendre(PARAM.whichmode,cos(pi*PARAM.remeshMapping(PARAM.t)));
P = P(1,:)';

%compute normal vector
[nx,ny] = normalVectorCurvilinear(x,y,PARAM.D1);

%perturb
x = x + nx.*P*delta;
y = y + ny.*P*delta;

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