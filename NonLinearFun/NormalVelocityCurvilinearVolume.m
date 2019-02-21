%compute velocity normal to the interface

function V = NormalVelocityCurvilinearVolume(dn,xBase,yBase,PARAM)

    %differentiation matrices
    D1 = PARAM.D1;

    %compute first derivative in the curvilinear parameter
    xpBase = D1*xBase;
    ypBase = D1*yBase;

    %normal vector
    [nxBase,nyBase] = normalVectorSplineSymmetric(xpBase,ypBase);
    
    %by symmetry
    nxBase([1 end]) = [1 -1];
    nyBase([1 end]) = [0 0];

    %cartesian coordiantes
    x = xBase+dn.*nxBase;
    y = yBase+dn.*nyBase;
    
    %compute solution
    V = axis_int_gauss_vect(x',y');

end