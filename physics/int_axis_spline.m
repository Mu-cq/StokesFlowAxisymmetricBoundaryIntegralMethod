%compute operator needed to integrate a quantity on a curved interface, the
%quantity in approximated linear on the elements and interface with cubic
%splines

function INT = int_axis_spline(x,y)

%compute spline
[~, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x,y);

%gauss points an weigths
GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152]';
GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170]';

%change metrics
GP = (GP+1)/2;
GW = GW/2;

%number of elements
I = numel(x)-1;
    
%many points and weigths
manyGP = repmat(GP,1,I);
manyGW = repmat(GW,1,I);
    
%many spline coordinates
bbx = repmat(bx,6,1);
ccx = repmat(cx,6,1);
ddx = repmat(dx,6,1);
aay = repmat(ay,6,1);
bby = repmat(by,6,1);
ccy = repmat(cy,6,1);
ddy = repmat(dy,6,1);
    
%for metrics in the integration
beta = aay+bby.*manyGP+ccy.*manyGP.^2+ddy.*manyGP.^3;
deta = bbx+2*ccx.*manyGP+3*ddx.*manyGP.^2;
dbeta = bby+2*ccy.*manyGP+3*ddy.*manyGP.^2;
h = sqrt(deta.^2+dbeta.^2);

%compute f in every gauss point assuming linear variation
phiA = 1-GP;    phiB = GP;
PHIA = repmat(phiA,1,I);
PHIB = repmat(phiB,1,I);
%manyPHI = fff(1:end-6).*PHIA + fff(7:end).*PHIB;

%compute operator for computing the integral
INT1 = sum(2*pi*PHIA.*beta.*h.*manyGW); %part for phiA
INT2 = sum(2*pi*PHIB.*beta.*h.*manyGW); %part for phiB
INT = [INT1 0] + [0 INT2];

end