

function varNew = volCorrBubbleBisection(t,var,Vin)

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
computeV = @(dn) axis_int_gauss_vect(xNew(dn)',yNew(dn)');

V0 = computeV(0);
R0 = nthroot(3/4*V0/pi,3);
n2 = 0.1*R0;    n1 = -0.1*R0;
errV = 1;
count = 1;
%bisection method
while abs(errV)>1e-14
    
    %compute displacement
    dn = (n1+n2)/2;
    
    %compute current volume
    V0 = computeV(dn);
    errV = (V0-Vin)/Vin;
    if errV>=0
        n2 = dn;
    elseif errV<0
        n1 = dn;
    end
       
    if count>1e3
        error('too many iterations')
    end
    count = count+1;
        
end

%compute neceassary normal displacement
%options = optimoptions('fsolve','TolFun',1e-14,'TolX',1e-14,'Display','none');
%dn = fsolve(computeV,0,options);

varNew(1:2:2*nNodes-1) = varNew(1:2:2*nNodes-1)+dn*nx;
varNew(2:2:2*nNodes) = varNew(2:2:2*nNodes)+dn*ny;