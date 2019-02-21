%compute vorticity from velocity field using finite difference 2nd oorder

function omega = vorticity(axial,radial,U,V,X,Y)

    %parametrix space
    dtx = 1/(axial-1);
    dty = 1/(radial-1);
    
    %diff matrix for vorticity
    DX = diag(ones(axial-1,1)/2,1) + diag(-ones(axial-1,1)/2,-1);    DX(1,1) = -1;   DX(1,2) = 1;   DX(end,end-1) = -1;   DX(end,end) = 1;
    DY = diag(ones(radial-1,1)/2,1) + diag(-ones(radial-1,1)/2,-1);    DY(1,1) = -1;    DY(1,2) = 1;    DY(end,end-1) = -1;    DY(end,end) = 1;
    
    %allocation
    dVdx = zeros(radial,axial);
    dUdy = zeros(radial,axial);
    
    %derivative in x direction
    for i = 1:radial
        
        %metrics terms
        x = X(i,:)';
        dXdt = DX*x/dtx;
        
        %derivative computation
        dVdt = DX*V(i,:)'/dtx;
        temp = dVdt./dXdt;
        dVdx(i,:) = temp';
        
    end
    
    %derivative in y direction
    for i = 1:axial
        
        %metrics terms
        y = Y(:,i);
        dYdt = DY*y/dty;
        
        %derivative computation
        dUdt = DY*U(:,i)/dtx;
        dUdy(:,i) = dUdt./dYdt;
        
    end
    
    %vorticity
    omega = dVdx-dUdy;

end