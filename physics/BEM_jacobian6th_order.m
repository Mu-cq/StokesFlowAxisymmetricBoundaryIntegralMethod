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
        aTemp5 = a; bTemp5 = b;
        aTemp5(i) = aTemp5(i) + 3*nrX(i)*dh;   bTemp5(i) = bTemp5(i) + 3*nrY(i)*dh;
        aTemp6 = a; bTemp6 = b;
        aTemp6(i) = aTemp6(i) - 3*nrX(i)*dh;   bTemp6(i) = bTemp6(i) - 3*nrY(i)*dh;
        
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
        
        [y,N] = bem_newton_extens(aTemp5,bTemp5,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u5 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        [y,N] = bem_newton_extens(aTemp6,bTemp6,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u6 = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (3/4*u1-3/4*u2-3/20*u3+3/20*u4+1/60*u5-1/60*u6)/dh;
        
    end

end