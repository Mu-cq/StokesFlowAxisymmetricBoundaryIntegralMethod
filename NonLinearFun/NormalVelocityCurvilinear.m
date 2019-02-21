%compute velocity normal to the interface

function [u,K1,K2,nx,ny] = NormalVelocityCurvilinear(dn,xBase,yBase,elem,lambda,capillary,PARAM)

    error('Correct bug in diff matrices')

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
    [sol,nx,ny,K1,K2] = bemNewtonExtensUp(x,y,elem,lambda,capillary,PARAM);
    
    %axial and radial
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    
    %compute normal velocity
    if PARAM.dropFrame==0
        u = nx.*ux + ny.*uy;
    elseif PARAM.dropFrame==1
        u = NormalVelocityDropFrame(x',y',ux,uy);
    end

end