%compute the volume integral of a solid of rotation around the axis

function V = axisIntGaussVectWeight(x,y,PARAM)

    %weights
    weightA = PARAM.weightA;
    weightB = PARAM.weightB;
    weightC = PARAM.weightC;
    weightD = PARAM.weightD;
    %weightSymA = PARAM.weightSymA;
    weightSymB = PARAM.weightSymB;
    weightSymC = PARAM.weightSymC;
    weightSymD = PARAM.weightSymD;
    
    xcm = center_mass(x,y);
    x = x-xcm;
    
    %ax = weightSymA*x';
    bx = weightSymB*x';
    cx = weightSymC*x';
    dx = weightSymD*x';
    ay = weightA*y';
    by = weightB*y';
    cy = weightC*y';
    dy = weightD*y';
    %[~, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
    
%     figure
%     plot(x,y,'o')
%     hold on
%     t = 0:1/1000:1;
%     for i = 1:numel(x)-1
%         xx = ax(i)+bx(i)*t+cx(i)*t.^2+dx(i)*t.^3;
%         yy = ay(i)+by(i)*t+cy(i)*t.^2+dy(i)*t.^3;
%         plot(xx,yy,'g-')
%     end
%     grid on
%     axis equal
%     hold off
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    GP = (GP+1)/2;
    GW = GW/2;
    
    I = numel(x)-1;
    
    manyGP = repmat(GP,1,I);
    manyGW = repmat(GW,1,I);
    
    bbx = reshape(repmat(bx,6,1),1,6*I);
    ccx = reshape(repmat(cx,6,1),1,6*I);
    ddx = reshape(repmat(dx,6,1),1,6*I);
    aay = reshape(repmat(ay,6,1),1,6*I);
    bby = reshape(repmat(by,6,1),1,6*I);
    ccy = reshape(repmat(cy,6,1),1,6*I);
    ddy = reshape(repmat(dy,6,1),1,6*I);
    
    beta = aay+bby.*manyGP+ccy.*manyGP.^2+ddy.*manyGP.^3;
    deta = bbx+2*ccx.*manyGP+3*ddx.*manyGP.^2;
        
    dV = -beta.*beta.*deta*pi*manyGW';
    
    V = sum(dV);

end