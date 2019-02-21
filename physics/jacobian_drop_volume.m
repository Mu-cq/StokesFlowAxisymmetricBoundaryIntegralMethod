%evaluate numerically the jacobian for the volume variation of a drop

function J = jacobian_drop_volume(a,b,V,dh)

    %initialize
    J = zeros(1,numel(a));
    
    %temp
    %r = sqrt(a.^2+b.^2);
    %q = numel(a)-1;
    %theta = 0:pi/q:pi;
    %INT = ([diff(theta),0]+[0,diff(theta)])/2;
    
    %vector in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(a)
        
        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(i) = aTemp(i) + nrX(i)*dh;   bTemp(i) = bTemp(i) + nrY(i)*dh;
        %rTemp = sqrt(aTemp.^2+bTemp.^2);
        
        Vnow = axis_int_gauss_vect(aTemp,bTemp);
        %Vnow = 2/3*pi*INT*(rTemp.^3.*sin(theta))';
        
        J(i) = (Vnow-V)/dh;
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end
    
end