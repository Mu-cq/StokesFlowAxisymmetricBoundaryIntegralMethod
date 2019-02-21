%gives the integration weights which perform a trapezi integral on the
%spline passing trough x and y

function weight = int_spline_symmetric_weight(x,y)

%compute spline
[~, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x,y);

%gauss points an weigths
GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];

%change metrics
GP = (GP+1)/2;
GW = GW/2;

%number of elements
I = numel(x)-1;
    
%many points and weigths
manyGP = repmat(GP,I,1);
manyGW = repmat(GW,I,1);
    
%many spline coordinates
bbx = repmat(bx,6,1)';
ccx = repmat(cx,6,1)';
ddx = repmat(dx,6,1)';
bby = repmat(by,6,1)';
ccy = repmat(cy,6,1)';
ddy = repmat(dy,6,1)';
    
%for metrics in the integration
%beta = aay+bby.*manyGP+ccy.*manyGP.^2+ddy.*manyGP.^3;
deta = bbx+2*ccx.*manyGP+3*ddx.*manyGP.^2;
dbeta = bby+2*ccy.*manyGP+3*ddy.*manyGP.^2;
h = sqrt(deta.^2+dbeta.^2);

%compute f in every gauss point assuming linear variation
phiA = 1-GP;    phiB = GP;

PHIA = repmat(phiA,I,1);
PHIB = repmat(phiB,I,1);

first = PHIA.*h.*manyGW;
second = PHIB.*h.*manyGW;

first = sum(first,2);
second = sum(second,2);

%compute integral
weight = [first; 0] + [0; second];
weight = weight';

end