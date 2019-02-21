%JacobianCM

function J = JacobianCM(a,b,XCM,dh)

    %initialize
    J = zeros(1,numel(a));

    %vector in radial direction
    nrX = a./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(b)

        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(i) = aTemp(i) + nrX(i)*dh;   bTemp(i) = bTemp(i) + nrY(i)*dh;
        
        XCMnow = center_mass(aTemp,bTemp);

        J(i) = (XCMnow-XCM)/dh; 

    end

end