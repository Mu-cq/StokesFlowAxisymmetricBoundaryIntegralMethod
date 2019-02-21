%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = conserveSurfaceAreaSpectralXY(xy,V0,Area0,PARAM)

%cartesian coordinates
x = xy(1:2:end-1);
y = xy(2:2:end);

%area constraint
DisEquality = [];
Equality = [surfaceCurvilinearAxisSpectral(x,y,PARAM)-Area0 VolumeCurvilinearAxisSpectral(x,y,PARAM)-V0];