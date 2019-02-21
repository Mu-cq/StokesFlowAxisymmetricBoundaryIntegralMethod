%integrate a certain quantity on an interface described by spline with
%Gauss integration

function F = int_axis_spline_symmetric(x,y,f)

%ok only if f is column
[~,col] = size(f);
if col>1
    error('Column vector has to be provided');
end

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
manyGP = repmat(GP,1,I);
manyGW = repmat(GW,1,I);
    
%many spline coordinates
bbx = reshape(repmat(bx,6,1),1,6*I);
ccx = reshape(repmat(cx,6,1),1,6*I);
ddx = reshape(repmat(dx,6,1),1,6*I);
aay = reshape(repmat(ay,6,1),1,6*I);
bby = reshape(repmat(by,6,1),1,6*I);
ccy = reshape(repmat(cy,6,1),1,6*I);
ddy = reshape(repmat(dy,6,1),1,6*I);
    
%for metrics in the integration
beta = aay+bby.*manyGP+ccy.*manyGP.^2+ddy.*manyGP.^3;
deta = bbx+2*ccx.*manyGP+3*ddx.*manyGP.^2;
dbeta = bby+2*ccy.*manyGP+3*ddy.*manyGP.^2;
h = sqrt(deta.^2+dbeta.^2);

%compute f in every gauss point assuming linear variation
phiA = 1-GP;    phiB = GP;
fff = reshape(repmat(f',6,1),1,6*I+6);
PHIA = reshape(repmat(phiA,I,1)',1,I*6);
PHIB = reshape(repmat(phiB,I,1)',1,6*I);
manyf = fff(1:end-6).*PHIA + fff(7:end).*PHIB;

%compute integral
F = 2*pi*manyf.*beta.*h*manyGW';

end