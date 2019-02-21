%compute numerical Jacobian by displacing interface

function J = BEM_jacobianLegendre(a,b,q,lambda,capillary,U,dh)

    %initialize
    J = zeros(numel(a));
    
    %radius
    r = sqrt(a.^2+b.^2);
    
    %independent varaibles
    theta = atan(b./a);
    theta = theta + pi*(theta<0);

    for i = 1:numel(a)
        
        %compute legemdre polynomial for perturbing
        %legendre function
        PPP = legendre(i-1,cos(theta));
        %legendre polynomia
        P = PPP(1,:);
        
        rTemp = r + dh*P;
        
        %display interface node in radial direction
        aTemp = rTemp.*cos(theta); bTemp = rTemp.*sin(theta);
        
        plot([aTemp flip(aTemp)],[bTemp -flip(bTemp)])
        axis equal
        grid on
        
        [y,N] = bem_newton_extens(aTemp,bTemp,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        %u = sqrt((N(1,:)'.*ux).^2 + (N(2,:)'.*uy).^2);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (u-U)./(dh*P)';
        
%         hold on
%         plot(u)
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end

end