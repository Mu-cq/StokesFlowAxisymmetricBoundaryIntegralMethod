%compute velocity normal to the interface

function [u,K1,K2,nx,ny] = NormalVelocityCurvilinearBond(dn,xBase,yBase,elem,lambda,capillary,Bond,PARAM)

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
    x = xBase+dn.*nxBase';
    y = yBase+dn.*nyBase';
    
    %compute solution
    [sol,nx,ny,K1,K2] = bemNewtonExtensUpBond(x,y,elem,lambda,capillary,Bond,PARAM);
    %[sol,~,~,nnx,nny,nnn] = BEM_Stokes(x,y,PARAM);
    
    %axial and radial
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    
    %compute normal velocity
    if PARAM.dropFrame==0
        u = nx.*ux + ny.*uy;
    elseif PARAM.dropFrame==1
        u = NormalVelocityDropFrame(x',y',ux,uy);
    end

end