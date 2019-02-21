%compute numerical Jacobian by displacing interface

function J = BEM_jacobianXY(a,b,q,lambda,capillary,U,dh)

    %initialize
    J = zeros(numel(a),2*numel(a));
    
    %jacobian displacing x
    for i = 1:numel(a)
        
        %display interface node in x direction
        aTemp = a; bTemp = b;
        aTemp(i) = aTemp(i) + dh;   %displace only x
        
        [y,N] = bem_newton_extens(aTemp,bTemp,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i) = (u-U)/dh;
        
    end
    
    %jacobian displacing y
    for i = 1:numel(a)
        
        %display interface node in y direction
        aTemp = a; bTemp = b;
        bTemp(i) = bTemp(i) + dh;   %displace only y
        
        [y,N] = bem_newton_extens(aTemp,bTemp,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        
        J(:,i+numel(a)) = (u-U)/dh;
        
    end

end