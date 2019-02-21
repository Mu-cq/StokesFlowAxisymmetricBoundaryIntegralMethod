%compute numerical Jacobian by displacing interface

function J = BEM_jacobian2Ddrop(a,b,PARAM,U,dh)

    n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
    aDrop = a(n+2*m+j+2:end);   bDrop = b(n+2*m+j+2:end);

    %initialize
    J = zeros(numel(aDrop)-1);

    %decide the displacement
    %dh = 1e-5;
    
    %xcm = center_mass(a,b);
    
    %vecor in radial direction
    nrX = aDrop./sqrt(aDrop.^2+bDrop.^2);    nrY = bDrop./sqrt(aDrop.^2+bDrop.^2);
    %nrX = (a-xcm)./sqrt(a.^2+b.^2);    nrY = b./sqrt(a.^2+b.^2);

    for i = 1:numel(aDrop)-1
        
        %display interface node in radial direction
        aTemp = a; bTemp = b;
        aTemp(n+2*m+j+1+i) = aTemp(n+2*m+j+1+i) + nrX(i)*dh;   bTemp(n+2*m+j+1+i) = bTemp(n+2*m+j+1+i) + nrY(i)*dh;
        
        if i==1
            aTemp(end) = aTemp(end) + nrX(i)*dh;    bTemp(end) = bTemp(end) + nrY(i)*dh;
        end
        
        %V = axis_int_gauss_vect(aTemp,bTemp);
        %alpha = (4/3*pi/V)^(1/3);
        %aTemp = alpha*aTemp;    bTemp = alpha*bTemp;
        
%         plot([aTemp flip(aTemp)],[bTemp -flip(bTemp)])
%         axis equal
%         grid on
        
        [y,~,N] = BEM_2DStokes_channel(aTemp,bTemp,PARAM);
        ux = y(2*n+4*m+2*j+1:2:end-1);  uy = y(2*n+4*m+2*j+2:2:end);
        u = N(1,:)'.*ux + N(2,:)'.*uy;
        u = u-u(1);
        
        J(:,i) = (u-U)/dh;
        
%         hold on
%         plot(u)
        
        %display(['Jacobian ' num2str(i) ' of ' num2str(numel(a)) ' is computed'])
        
    end

end