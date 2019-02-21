%compute integration weights for integration on a axisymmetric line

function weight = integrationOnLineWeightAxis(x,y,orderVariable,orderGeometry,SPlinesType)

if orderVariable==0 && orderGeometry==0             % in the midlle of the straight element
   
    dl = sqrt(diff(x).^2+diff(y).^2);
    Ymiddle = (y(1:end-1)+y(2:end))/2;
    weight = 2*pi*dl.*Ymiddle;
    
elseif orderVariable==1 && orderGeometry==0         % straight element
   
    dl = sqrt(diff(x).^2+diff(y).^2);
    weight = pi*([y(1:end-1).*dl 0]+[0 y(2:end).*dl]);
    
elseif orderVariable==0 && orderGeometry==1     % in the midlle of the curved element
    
    %compute spline coeff
    if SPlinesType==1
         [~,bx,cx,dx,ay,by,cy,dy] = spline_natural(x,y);
    elseif SPlinesType==2
         [~,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x,y);
    end
    
    %integration weights
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %adapt for interval between [0,1]
    GP = (GP+1)/2;
    GW = GW/2;
    
    %prepare coeff
    manyGP = repmat(GP,numel(bx),1);
    manyGW = repmat(GW,numel(bx),1);
    bx = repmat(bx',1,6);
    cx = repmat(cx',1,6);
    dx = repmat(dx',1,6);
    ay = repmat(ay',1,6);
    by = repmat(by',1,6);
    cy = repmat(cy',1,6);
    dy = repmat(dy',1,6);
    
    %compute arc lenght
    xp = derSplinesUnk(manyGP,bx,cx,dx);
    yp = derSplinesUnk(manyGP,by,cy,dy);
    yAll = valSplinesUnk(manyGP,ay,by,cy,dy);
    I = manyGW.*yAll.*sqrt(xp.^2+yp.^2);
    weight = 2*pi.*sum(I,2)';
    
elseif orderVariable==1 && orderGeometry==1     % on the sides of the curved element
    
    %compute spline coeff
    if SPlinesType==1
         [~,bx,cx,dx,ay,by,cy,dy] = spline_natural(x,y);
    elseif SPlinesType==2
         [~,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x,y);
    end
    
    %integration weights
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %adapt for interval between [0,1]
    GP = (GP+1)/2;
    GW = GW/2;
    
    %linear test function
    phia = 1-GP;
    phib = GP;
    PHIA = repmat(phia,numel(bx),1);
    PHIB = repmat(phib,numel(bx),1);
    
    %prepare coeff
    manyGP = repmat(GP,numel(bx),1);
    manyGW = repmat(GW,numel(bx),1);
    bx = repmat(bx',1,6);
    cx = repmat(cx',1,6);
    dx = repmat(dx',1,6);
    ay = repmat(ay',1,6);
    by = repmat(by',1,6);
    cy = repmat(cy',1,6);
    dy = repmat(dy',1,6);
    
    %compute arc lenght
    xp = derSplinesUnk(manyGP,bx,cx,dx);
    yp = derSplinesUnk(manyGP,by,cy,dy);
    yAll = valSplinesUnk(manyGP,ay,by,cy,dy);
    IA = manyGW.*yAll.*sqrt(xp.^2+yp.^2).*PHIA;
    IB = manyGW.*yAll.*sqrt(xp.^2+yp.^2).*PHIB;
    weightA = sum(IA,2)';
    weightB = sum(IB,2)';
    weight = 2*pi*([weightA 0]+[0 weightB]);
    
end













