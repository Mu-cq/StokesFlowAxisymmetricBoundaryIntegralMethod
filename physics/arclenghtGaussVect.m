%compute surface integral with trapezi rule in cylindrical coordiante

function [l,dl] = arclenghtGaussVect(x,y)

    %using trapezi rules for integration
    %dtheta = pi/(numel(x)-1);
    
%     a = find(x==0);
%     b = find(y==0);
%     
%     x = x(1:a(1)-1);
%     y = y(1:b(2));
    
    xcm = center_mass(x,y);
    x = x-xcm;
    
    [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(x, y);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152]';
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170]';
    
    GP = (GP+1)/2;
    GW = GW/2;

    I = numel(x)-1;
    
    manyGP = repmat(GP,1,I);
    manyGW = repmat(GW,1,I);
    
    bbx = repmat(bx,6,1);
    ccx = repmat(cx,6,1);
    ddx = repmat(dx,6,1);
    bby = repmat(by,6,1);
    ccy = repmat(cy,6,1);
    ddy = repmat(dy,6,1);
    
    deta = bbx+2*ccx.*manyGP+3*ddx.*manyGP.^2;
    dbeta = bby+2*ccy.*manyGP+3*ddy.*manyGP.^2;
        
    dl = sum(sqrt(deta.*deta+dbeta.*dbeta).*manyGW);
    
    l = [0 cumsum(dl)];

end