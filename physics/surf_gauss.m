%compute surface integral with trapezi rule in cylindrical coordiante

function [A,dA] = surf_gauss(x,y)

    %using trapezi rules for integration
    %dtheta = pi/(numel(x)-1);
    
%     a = find(x==0);
%     b = find(y==0);
%     
%     x = x(1:a(1)-1);
%     y = y(1:b(2));
    
    xcm = center_mass(x,y);
    x = x-xcm;
    
    [~, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    GP = (GP+1)/2;
    GW = GW/2;

    I = numel(x)-1;
    
    dA = zeros(I,1);
    
    for i = 1:I
        %eta = ax(i)+bx(i)*GP(k)+cx(i)*GP(k)^2+dx(i)*GP(k)^3;
        beta = ay(i)+by(i)*GP+cy(i)*GP.^2+dy(i)*GP.^3;
        deta = bx(i)+2*cx(i)*GP+3*dx(i)*GP.^2;
        dbeta = by(i)+2*cy(i)*GP+3*dy(i)*GP.^2;
        
        dA(i) = (beta.*sqrt(deta.*deta+dbeta.*dbeta)*GW')*2*pi;
    end
    
    A = sum(dA);

end