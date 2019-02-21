%compute numerical Jacobian by displacing interface

function J2 = BEM_jacobian2nd_derivative_diss(a,b,q,lambda,capillary,dh,J1)

    %initialize
    %J2 = zeros(numel(a)^2,numel(a));
    J2 = zeros(numel(a),numel(a),numel(a));
    
    %because is econd order
    %dh = 2*dh;
    
    %vecor in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(a)
        
        %displace interface node in radial direction
        aTemp1 = a; bTemp1 = b;
        aTemp1(i) = aTemp1(i) + nrX(i)*dh;   bTemp1(i) = bTemp1(i) + nrY(i)*dh;
        
        %compute Jacobian of the perturbed interface
        [y,~,~,~,~,~,~,~,~,~,FX,FY] = bem_newton_extens(aTemp1,bTemp1,q,lambda,capillary);
        ux = y(1:2:end-1);  uy = y(2:2:end);
        G = FX'.*ux + FY'.*uy;
        J = BEM_jacobian(aTemp1,bTemp1,q,lambda,capillary,G,dh);
        
        %compute jacobian of the jacobian (second derivative of the velocity)
        %J2(1+(i-1)*numel(a):i*numel(a),:) = (J-J1)/dh;
        J2(:,:,i) = (J-J1)/dh;
        
    end
    
    %J2 = reshape(J2,numel(a),numel(a)^2);

end