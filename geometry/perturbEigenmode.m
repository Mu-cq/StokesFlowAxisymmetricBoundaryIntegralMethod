%pertub shape with second legendre polynomial

function [x,y] = perturbEigenmode(x,y,delta,PARAM)

display(['Perturb shape with ' num2str(PARAM.whichmode) 'th most unstable mode, amplitude=' num2str(delta)])

%from grid to modes
[xMode,yMode] = fromGridToModes(x,y,PARAM);
xMode = xMode(1:PARAM.dealiasing);
yMode = yMode(1:PARAM.dealiasing);

%current center of mass
xcm = CenterMassCurvAxisSpectral(x,y,PARAM);

%nonlinear function
if PARAM.BC==1
    fNonlinear = @(unk) NormalVelocitySpectralVolumeXY2modes(unk,xMode,yMode,PARAM.V0,xcm,PARAM);
elseif PARAM.BC==2
   error('Not implemented')
end

%perturbation
if PARAM.dropFrame==0
    perturb = zeros(numel(xMode)-1,1);
elseif PARAM.dropFrame==1 || PARAM.dropFrame==2
    perturb = zeros(numel(xMode)-2,1);
end

%Compute Jacobian
J = JacobianHandle(fNonlinear,perturb,PARAM.dh);

%compute eigenvalues and eigenfunctions
[eigenmode,eigenval] = eig(J);
eigenval  = diag(eigenval);

%order from small to big
[~,ind] = sort(real(eigenval),'descend');
eigenmode = eigenmode(:,ind);

%choose mode
mode = eigenmode(:,PARAM.whichmode);
%mode(PARAM.dealiasing+1:numel(x)) = 0;

%compute physical eigenmode
if PARAM.legendre==1||PARAM.legendre==2
     PPP = PARAM.PPP;
elseif PARAM.legendre==0
     PPP = PARAM.TTT;
end

%perturbation of the physical grid
if PARAM.dropFrame==0
        MODE = sum(repmat(mode,1,numel(x)).*PPP(2:PARAM.dealiasing,:))';
elseif PARAM.dropFrame==1 || PARAM.dropFrame==2
        MODE = sum(repmat(mode,1,numel(x)).*PPP(3:PARAM.dealiasing,:))';
end

%compute normal vector
[nx,ny] = normalVectorCurvilinear(x,y,PARAM.D1);

%perturb
x = x + nx.*MODE*delta/MODE(1);
y = y + ny.*MODE*delta/MODE(1);

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


