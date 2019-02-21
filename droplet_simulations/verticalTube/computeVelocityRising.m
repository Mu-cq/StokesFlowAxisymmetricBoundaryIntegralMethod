% compute velocity of the motor and bubble once the geomtry is known

function [Udrop,yStokes,x,y,PARAM,Vdrop] = computeVelocityRising(t,var,tParametricBase,PARAM)

%bubble coordinates
nNodesBubble = numel(var)/2;
xDrop = var(1:2:2*nNodesBubble-1);
yDrop = var(2:2:2*nNodesBubble);
yDrop([1 end]) = [0 0];

%compute curvilinear bubble
dl = sqrt(diff(xDrop).^2+diff(yDrop).^2);
l = [0 cumsum(dl')];
PARAM.ppx(4) = spline(l/l(end),xDrop');
PARAM.ppy(4) = spline(l/l(end),yDrop');

%remesh
[xBase,yBase] = buildGeometryPanelsParametric(tParametricBase,PARAM);
xBase = xBase(1:3); xBase{4} = xDrop';
yBase = yBase(1:3); yBase{4} = yDrop';

tParametric = remeshPanels(xBase,yBase,tParametricBase,tParametricBase,PARAM);
for k = 1:numel(PARAM.n)
    PARAM.n(k) = numel(tParametric{k})-1;
end

%build shape
[xBlock1,yBlock1] = buildGeometryPanelsParametric(tParametric,PARAM);
xBlock1 = xBlock1(1:3);
yBlock1 = yBlock1(1:3);
x = xBlock1; x{4} = xDrop';
y = yBlock1; y{4} = yDrop';
PARAM.n(4) = nNodesBubble-1;

figure(100)
hold off
%plotGeometryDrop(x,y,PARAM,1)
plotGeometryStokes(x,y,0,[],[],[],0,PARAM)
plot(x{4},y{4},'xr')
plot(x{2},y{2},'xr')
axis([-4 4 -1 1])
xlabel('x')
ylabel('r')
drawnow

%solve Stokes equation
[yStokes,~,~,nnx,nny,nnn] = BEM_Stokes(x,y,PARAM);

%drop velocities
Udrop = yStokes(2*sum(nnn(1:3))+1:end);
Ux = Udrop(1:2:end-1);
Uy = Udrop(2:2:end);

%drop velocity
Vdrop = DropVelocityAxis(x{4},y{4},Ux.*nnx{4}'+Uy.*nny{4}');

%project in the normal direction
Udrop = zeros(2*numel(Ux),1);
Un = (Ux-Vdrop).*nnx{4}' + Uy.*nny{4}';
Udrop(1:2:end-1) = Un.*nnx{4}';
if PARAM.dropFrame==0
    Udrop = Udrop+Vdrop;
end
Udrop(2:2:end) = Un.*nny{4}';







