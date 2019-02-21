%compute velocity normal to the interface

function V = volumeCurvilinear2(dn,xBase,yBase)

    %compute first derivative in the curvilinear parameter
    [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(xBase',yBase');
    xpBase = derSplines(bx,cx,dx)';
    ypBase = derSplines(by,cy,dy)';

    %normal vector
    [nxBase,nyBase] = normalVectorSplineSymmetric(xpBase,ypBase);
    
    %by symmetry
    nxBase([1 end]) = [1 -1];
    nyBase([1 end]) = [0 0];

    %cartesian coordnantes
    x = xBase+dn.*nxBase;
    y = yBase+dn.*nyBase;
    
    %compute solution
    V = axis_int_gauss_vect(x',y');

end