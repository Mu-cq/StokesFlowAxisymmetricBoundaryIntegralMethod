%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = conserveSurfaceAreaSpectral(perturb,xBase,yBase,V0,Area0,PARAM)

%cartesian coordinates
[x,y] = getVolumeConservingShapeSpectral(perturb,xBase,yBase,V0,PARAM);

%area constraint
DisEquality = [];
Equality = surfaceCurvilinearAxisSpectral(x,y,PARAM)-Area0;