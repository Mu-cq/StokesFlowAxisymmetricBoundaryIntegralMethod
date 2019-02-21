%compute normal vector using symmetric splines

function [nx,ny] = normalVectorSplineSymmetric(xp,yp)

h = sqrt(xp.^2+yp.^2);

nx = yp./h;
ny = -xp./h;