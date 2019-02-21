%compute numerical Jacobian by displacing interface

function J = BEM_jacobian(a,b,q,lambda,capillary,U,dh)

    %initialize
    J = zeros(numel(a));

    %decide the displacement
    %dh = 1e-5;
    
    %xcm = center_mass(a,b);
    
    %vecor in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);
    %nrX = (a-xcm)./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(a)
        
        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(i) = aTemp(i) + nrX(i)*dh;   bTemp(i) = bTemp(i) + nrY(i)*dh;
        
%         V = axis_int_gauss_vect(aTemp,bTemp);
%         alpha = (4/3*pi/V)^(1/3);
%         aTemp = alpha*aTemp;    bTemp = alpha*bTemp;
        
%         plot([aTemp flip(aTemp)],[bTemp -flip(bTemp)])
%         axis equal
%         grid on
        
        [y,N] = bem_newton_extens(aTemp,bTemp,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        %u = sqrt((N(1,:)'.*ux).^2 + (N(2,:)'.*uy).^2);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (u-U)/dh;
        
%         hold on
%         plot(u)
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end

end