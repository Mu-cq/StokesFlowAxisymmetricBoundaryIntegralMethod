%compute velocity normal to the interface

function xcm = centerOfMassCurvilinear2(dn,xBase,yBase)

    %compute spline coeff
    [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(xBase',yBase');
    
    %first and second derivative with splines
    xpBase = derSplines(bx,cx,dx);
    ypBase = derSplines(by,cy,dy);

    %normal vector
    [nxBase,nyBase] = normalVectorSplineSymmetric(xpBase,ypBase);
    
    %by symmetry
    nxBase([1 end]) = [1 -1];
    nyBase([1 end]) = [0 0];

    %cartesian coordiantes
    x = xBase+dn.*nxBase;
    y = yBase+dn.*nyBase;
    
    %compute solution
    xcm = center_mass(x',y');

end