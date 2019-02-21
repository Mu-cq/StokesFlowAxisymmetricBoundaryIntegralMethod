%aadd the motion of the center of mass

function xyMode = addMotionCenterMass(dt,xyMode,Vdrop,PARAM)

xMode = xyMode(1:2:end-1);
yMode = xyMode(2:2:end);

%from modes to nodal value
[x,y] = fromModesToGrid(xMode,yMode,PARAM);

%displace
x = x + dt*Vdrop;

%compute modes
[xMode,yMode] = fromGridToModes(x,y,PARAM);

xyMode(1:2:end-1) = xMode;
xyMode(2:2:end) = yMode;