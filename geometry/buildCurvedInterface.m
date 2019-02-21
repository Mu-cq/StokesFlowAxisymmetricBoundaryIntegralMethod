%build curved interface

function [x,y] = buildCurvedInterface(xDot,yDot,nElem)

%fit with spline
xTemp = linspace(xDot(1),xDot(end),1e3);
yTemp = spline(xDot,yDot,xTemp);

%compute arc lenght
dx = diff(xTemp);
dy = diff(yTemp);
dl = sqrt(dx.^2+dy.^2);
l = [0 cumsum(dl)];

x = spline(l,xTemp,linspace(0,l(end),nElem+1));
y = spline(l,yTemp,linspace(0,l(end),nElem+1));

%smooth
[~,x] = spaps(linspace(0,l(end),nElem+1),x,1);
[~,y] = spaps(linspace(0,l(end),nElem+1),y,1);

figure
plot(x,y,'-x')
hold on
plot(xDot,yDot,'o')
xlabel('x')
ylabel('y')
grid on
axis equal
drawnow