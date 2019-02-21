% compute velocity of the motor and bubble once the geomtry is known

function [Udrop,yStokes,x,y,PARAM,Vdrop] = computeVelocityDrop(t,var,PARAM)

%bubble coordinates
nNodesBubble = numel(var)/2;
xDrop = var(1:2:2*nNodesBubble-1);
yDrop = var(2:2:2*nNodesBubble);
yDrop([1 end]) = [0 0];

%build shape
x{1} = xDrop';
y{1} = yDrop';
PARAM.n(1) = nNodesBubble-1;

% figure(100)
% hold off
% %plotGeometryDrop(x,y,PARAM,1)
% plotGeometryStokes(x,y,0,[],[],[],0,PARAM)
% plot(x{1},y{1},'xr')
% axis([-4 4 -1.5 1.5])
% xlabel('x')
% ylabel('r')
% drawnow

%solve Stokes equation
[yStokes,~,~,nnx,nny] = BEM_Stokes(x,y,PARAM);

%drop velocities
Udrop = yStokes;
Ux = Udrop(1:2:end-1);
Uy = Udrop(2:2:end);

%drop velocity
Vdrop = DropVelocityAxis(x{1},y{1},Ux.*nnx{1}'+Uy.*nny{1}');

%project in the normal direction
Udrop = zeros(2*numel(Ux),1);
Un = (Ux-Vdrop).*nnx{1}' + Uy.*nny{1}';
Udrop(1:2:end-1) = Un.*nnx{1}';
if PARAM.dropFrame==0
    Udrop = Udrop+Vdrop;
end
Udrop(2:2:end) = Un.*nny{1}';







