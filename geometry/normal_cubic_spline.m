function [hm,nx,ny ] = normal_cubic_spline( x,y )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nn = length(x); %number of nodes
ne = nn - 1;    %number of elements
K  = ne - 1;    %size of linear system  eq 3.1.27                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        lenelem = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);
lennode = x;
lennode(1) = 0;

a = x(1:K); s = a;
b = a(1:(end-1));
c = b;

lenelem = sqrt((x(2:end)-x(1:(end-1))).^2+(y(2:end)-y(1:(end-1))).^2);
lenelhf = 0.5*lenelem;

for inode = 2:nn
    lennode(inode) = lennode(inode-1)+lenelem(inode-1);
end

bx = x;
bx(1)=0;bx(end)=0;%%boundary condition for natural cubic spline.
by = bx;

a = 2./3*(lenelem(1:K)+lenelem(2:(K+1)));
b = 1./3*lenelem(2:K);
c = b;

sx = (x(3:(K+2))-x(2:(K+1)))./lenelem(2:(K+1)) - (x(2:(K+1))-x(1:K))./lenelem(1:K);
sy = (y(3:(K+2))-y(2:(K+1)))./lenelem(2:(K+1)) - (y(2:(K+1))-y(1:K))./lenelem(1:K);
[ xtmp ] = thomas_linear( a,b,c,sx );
bx(2:(end-1)) = xtmp;
[ xtmp ] = thomas_linear( a,b,c,sy );
by(2:(end-1)) = xtmp;

[ ax,cx ] = update_ac_b(bx,lenelem,x,ne);
[ ay,cy ] = update_ac_b(by,lenelem,y,ne);


nx = (3*ay.*lenelhf.^2+2*by(1:(end-1)).*lenelhf+cy);
ny = (3*ax.*lenelhf.^2+2*bx(1:(end-1)).*lenelhf+cx);
hm = sqrt(nx.^2+ny.^2);
nx = -1./hm.*nx;
ny =  1./hm.*ny;
end

