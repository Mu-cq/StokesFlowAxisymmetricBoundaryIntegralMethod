%Stokeslet 3D

function [GXX,GXY,GXZ,GYX,GYY,GYZ,GZX,GZY,GZZ] = axisKernelSingleLayer3D(x,y,z,x0,y0,z0)

dx = x0-x;
dy = y0-y;
dz = z0-z;
r = sqrt(dx.^2+dy.^2+dz.^2);

%Stokeslet
GXX = 1./r + dx.^2./r.^3;
GXY = dx.*dy./r.^3;
GXZ = dx.*dz./r.^3;
GYX = GXY;
GYY = 1./r + dy.^2./r.^3;
GYZ = dy.*dz./r.^3;
GZX = GXZ;
GZY = GYZ;
GZZ = 1./r + dz.^2./r.^3;