%compute normal vector

function dl = computeDL(x,y,orderVariable,orderGeometry,SPlinesType)

if (orderVariable==0||orderVariable==1) && orderGeometry==0         % in the midlle of the straight element
   
    dl = sqrt(diff(x).^2+diff(y).^2);
    
elseif (orderVariable==0||orderVariable==1) && orderGeometry==1     % in the midlle of the curved element
    
    %compute spline coeff
    if SPlinesType==1
         [~,bx,cx,dx,~,by,cy,dy] = spline_natural(x,y);
    elseif SPlinesType==2
         [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(x,y);
    end
    
    %first and second derivative with splines
    xp = @(t) derSplinesUnk(t,bx,cx,dx);
    yp = @(t) derSplinesUnk(t,by,cy,dy);
    
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
    by = repmat(by',1,6);
    cy = repmat(cy',1,6);
    dy = repmat(dy',1,6);
    
    %compute arc lenght
    xp = derSplinesUnk(manyGP,bx,cx,dx);
    yp = derSplinesUnk(manyGP,by,cy,dy);
    I = manyGW.*sqrt(xp.^2+yp.^2);
    dl = sum(I,2)';
    
end













