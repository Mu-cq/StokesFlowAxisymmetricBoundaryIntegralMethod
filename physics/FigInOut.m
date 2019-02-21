%figure if the points are inside or outside the droplet with a line
%integral: if phi=2*pi is inside, if phi=0 is outside, if phi=pi is on the
%interface. In an approximate way I can hypotheze that it will never be on
%the interface

function phi = FigInOut(a,b,X0,Y0)

    el = numel(a)-1;
    
    x = [a a(end-1:-1:2)];
    y = [b -b(end-1:-1:2)];

    X0matr = repmat(X0,1,2*el-1);
    Y0matr = repmat(Y0,1,2*el-1);
    
    globalA = repmat(x,numel(X0),1);
    globalB = repmat(y,numel(X0),1);
    
    dx1 = globalA(:,1:end-1)-X0matr;
    dy1 = globalB(:,1:end-1)-Y0matr;
    dx2 = globalA(:,2:end)-X0matr;
    dy2 = globalB(:,2:end)-Y0matr;
    
    rr = sqrt(dx1.*dx1+dy1.*dy1).*sqrt(dx2.*dx2+dy2.*dy2);

    A = asin((dy2.*dx1-dx2.*dy1)./rr);
    
    phi = sum(A,2);
    
end