%compute numerical Jacobian by displacing interface

function J = BEM_jacobian_dissipation(a,b,q,lambda,capillary,G1,dh)

    %initialize
    J = zeros(numel(a));
    
    %vecor in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(a)
        
        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(i) = aTemp(i) + nrX(i)*dh;   bTemp(i) = bTemp(i) + nrY(i)*dh;
        
        [y,~,~,~,~,~,~,~,~,~,fx,fy] = bem_newton_extens(aTemp,bTemp,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        G = ux.*fx' + uy.*fy';
        
        J(:,i) = (G-G1)/dh;
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end

end