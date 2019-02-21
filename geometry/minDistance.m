%compute minimum distnce betweeen points

function dl = minDistance(x,y,z,x0,y0,z0)

N = numel(x);
M = numel(x0);

xx =repmat(x,1,M);
yy =repmat(y,1,M);
zz =repmat(z,1,M);
xx0 =repmat(x0',N,1);
yy0 =repmat(y0',N,1);
zz0 =repmat(z0',N,1);

%distance
dx2 = (xx-xx0).^2;
dy2 = (yy-yy0).^2;
dz2 = (zz-zz0).^2;
dl = sqrt(dx2+dy2+dz2);

%minimum distance
for i = 1:N
    for k = 1:M
        
        if dl(i,k)==0
            dl(i,k) = 10;
        end
        
    end
end
dl = min(min(dl));