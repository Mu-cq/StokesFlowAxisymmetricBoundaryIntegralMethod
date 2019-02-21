%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = conserveSurfaceAreaBEM(XY,V0,Area0)

%cartesian coordinates
x = XY(1:2:end-1);
y = XY(2:2:end);

%area constraint
DisEquality = [];
Equality = [surf_gauss_vect(x',y')-Area0 axis_int_gauss_vect(x',y')-V0];