%compute numerical Jacobian by displacing interface

function J = BEM_jacobian(a,b,q,lambda,capillary,U,dh)

    %initialize
    J = zeros(numel(a));

    %vecor in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(a)
        
        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(i) = aTemp(i) + nrX(i)*dh;   bTemp(i) = bTemp(i) + nrY(i)*dh;
        
        [y,N] = bem_newton_extens(aTemp,bTemp,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (u-U)/dh;
        
    end

end