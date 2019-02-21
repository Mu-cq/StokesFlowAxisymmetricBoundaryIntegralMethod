%compute numerical Jacobian by displacing interface

function J = BEM_jacobian4th_order(a,b,q,lambda,capillary,dh)

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
        aTemp3 = a; bTemp3 = b;
        aTemp3(i) = aTemp3(i) + 2*nrX(i)*dh;   bTemp3(i) = bTemp3(i) + 2*nrY(i)*dh;
        aTemp4 = a; bTemp4 = b;
        aTemp4(i) = aTemp4(i) - 2*nrX(i)*dh;   bTemp4(i) = bTemp4(i) - 2*nrY(i)*dh;
        
        [y,N] = bem_newton_extens(aTemp1,bTemp1,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u1 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        [y,N] = bem_newton_extens(aTemp2,bTemp2,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u2 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        [y,N] = bem_newton_extens(aTemp3,bTemp3,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u3 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        [y,N] = bem_newton_extens(aTemp4,bTemp4,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u4 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (2/3*u1-2/3*u2-1/12*u3+1/12*u4)/dh;
        
    end

end