%compute numerical Jacobian by displacing interface

function J = BEM_jacobianDropChannel(a,b,PARAM,U,dh)

    n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
    aDrop = a(n+m+j+2:end);   bDrop = b(n+m+j+2:end);

    %initialize
    J = zeros(numel(aDrop));

    %decide the displacement
    %dh = 1e-5;
    
    %xcm = center_mass(a,b);
    
    %vecor in radial direction
    nrX = aDrop./sqrt(aDrop.^2+bDrop.^2);    nrY = bDrop./sqrt(aDrop.^2+bDrop.^2);
    %nrX = (a-xcm)./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(aDrop)
        
        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(n+m+j+1+i) = aTemp(n+m+j+1+i) + nrX(i)*dh;   bTemp(n+m+j+1+i) = bTemp(n+m+j+1+i) + nrY(i)*dh;
        
        %V = axis_int_gauss_vect(aTemp,bTemp);
        %alpha = (4/3*pi/V)^(1/3);
        %aTemp = alpha*aTemp;    bTemp = alpha*bTemp;
        
%         plot([aTemp flip(aTemp)],[bTemp -flip(bTemp)])
%         axis equal
%         grid on
        
        [y,~,~,~,N] = BEM_drop_channel(aTemp,bTemp,PARAM);
        ux = y(2*n+2*m+2*j+1:2:end-1);  uy = y(2*n+2*m+2*j+2:2:end);
        u = N(1,:)'.*(ux-ux(1)) + N(2,:)'.*uy;
        
        J(:,i) = (u-U)/dh;
        
%         hold on
%         plot(u)
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end

end