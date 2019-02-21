%compute numerical Jacobian by displacing interface

function J = BEM_jacobianManyAlpha2(theta,r,alpha,q,lambda,capillary,U,dh)

    %initialize
    J = zeros(numel(r));

    for i = 1:numel(r)
        
        %display interface node in radial direction
        rTemp = r;
        alpha(i) = alpha(i) + dh;
        a = rTemp'.*cos(theta).*alpha';     b = rTemp'.*sin(theta).*alpha';
                
        [y,N] = bem_newton_extens(a,b,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (u-U)/dh;
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end

end