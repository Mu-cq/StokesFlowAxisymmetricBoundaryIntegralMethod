%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = conserveSurfaceAreaSpectralXYmodes(xy,V0,Area0,PARAM)

xMode = xy(1:2:end-1);
yMode = xy(2:2:end);

%cartesian coordinates
[x,y] = fromModesToGrid(xMode,yMode,PARAM);

%area constraint
DisEquality = [];
Equality = [surfaceCurvilinearAxisSpectral(x,y,PARAM)-Area0 VolumeCurvilinearAxisSpectral(x,y,PARAM)-V0];