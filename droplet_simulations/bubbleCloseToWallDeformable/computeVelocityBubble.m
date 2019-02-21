% compute velocity of the motor and bubble once the geomtry is known

function [Udrop,yStokes,x,y,PARAM] = computeVelocityBubble(t,var,PARAM)

%bubble coordinates
nNodesBubble = numel(var)/2;
xDrop = var(1:2:2*nNodesBubble-1);
yDrop = var(2:2:2*nNodesBubble);
yDrop([1 end]) = [0 0];
x{1} = xDrop';
y{1} = yDrop';
PARAM.n = numel(xDrop)-1;

figure(100)
%hold off
plot(xDrop,yDrop,'k')
hold on
plot(xDrop,-yDrop,'k','MarkerSize',25)
plot([0 0],[-10 10],'k')
axis equal
axis([0 2*max(xDrop) -2*max(yDrop) 2*max(yDrop)])
xlabel('$z$','interpreter','latex')
ylabel('$r$','interpreter','latex')
%camroll(90)
drawnow

%solve Stokes equation
[yStokes,~,~,nnx,nny] = BEM_Stokes(x,y,PARAM);

%drop velocities
Udrop = yStokes;
Ux = Udrop(1:2:end-1);
Uy = Udrop(2:2:end);

%use only normal velocity
if PARAM.Unormal==1
    
    Un = Ux.*nnx{1}' + Uy.*nny{1}';
    Ux = Un.*nnx{1}';
    Uy = Un.*nny{1}';
    Udrop(1:2:end-1) = Ux;
    Udrop(2:2:end) = Uy;
    
end







