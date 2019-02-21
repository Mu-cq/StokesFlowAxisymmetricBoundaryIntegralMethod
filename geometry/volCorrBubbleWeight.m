

function varNew = volCorrBubbleWeight(t,var,Vin,PARAM)

%initialize
varNew = var;

%bubble shape
nNodes = floor(numel(var)/2);
xBubble = var(1:2:2*nNodes-1);
yBubble = var(2:2:2*nNodes);

%compute splines coeff
[~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(xBubble',yBubble');
    
%compute normal vector
[nx,ny] = DropNormalVector(bx,cx,dx,by,cy,dy);

%compute new coordinates
xNew = @(dn) xBubble+dn*nx;
yNew = @(dn) yBubble+dn*ny;

%compute volume
computeV = @(dn) (axisIntGaussVectWeight(xNew(dn)',yNew(dn)',PARAM)-Vin)/Vin;

%compute neceassary normal displacement
options = optimoptions('fsolve','TolFun',1e-14,'TolX',1e-14,'Display','none');
dn = fsolve(computeV,0,options);

varNew(1:2:2*nNodes-1) = varNew(1:2:2*nNodes-1)+dn*nx;
varNew(2:2:2*nNodes) = varNew(2:2:2*nNodes)+dn*ny;