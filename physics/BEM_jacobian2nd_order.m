%compute numerical Jacobian by displacing interface

function J = BEM_jacobian2nd_order(a,b,q,lambda,capillary,dh)

    %initialize
    J = zeros(numel(a));
    
    %vecor in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(a)
        
        %display interface node in radial direction
        aTemp1 = a; bTemp1 = b;
        aTemp1(i) = aTemp1(i) + nrX(i)*dh;   bTemp1(i) = bTemp1(i) + nrY(i)*dh;
        aTemp2 = a; bTemp2 = b;
        aTemp2(i) = aTemp2(i) - nrX(i)*dh;   bTemp2(i) = bTemp2(i) - nrY(i)*dh;
        
        [y,N] = bem_newton_extens(aTemp1,bTemp1,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u1 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        [y,N] = bem_newton_extens(aTemp2,bTemp2,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u2 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (0.5*u1-0.5*u2)/dh;
        
    end

end