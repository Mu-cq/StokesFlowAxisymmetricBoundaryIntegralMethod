%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_2Dstokes_linear(x,y,normal)

    %tic
    
    %number of singularities
    N = numel(x)-1;
    el = N;
    
    %compute normals
    DX = x(2:end)-x(1:end-1);
    DY = y(2:end)-y(1:end-1);
    
    if normal==1
        r = [DY./sqrt(DX.*DX+DY.*DY); -DX./sqrt(DX.*DX+DY.*DY)];
    elseif normal==0
        r = -[DY./sqrt(DX.*DX+DY.*DY); -DX./sqrt(DX.*DX+DY.*DY)];
    end
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x(1:end-1);
    Y0 = y(1:end-1);
    
    %moltiplicate sing points coordinates
    tempx = repmat(X0,6,1);
    XX = reshape(tempx,1,6*N);
    tempy = repmat(Y0,6,1);
    YY = reshape(tempy,1,6*N);
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %modify GP
    GP = GP+1;
    
    %every Gauss point
    GPX = repmat(GP,1,el).*reshape((repmat(deltaX/2,6,1)),1,6*el);
    GPY = repmat(GP,1,el).*reshape((repmat(deltaY/2,6,1)),1,6*el);
    
    phia = 1-GP/2;
    phib = GP/2;
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %compute green's function
    [SXX,SXY,SYY,QXXX,QXXY,QXYY,QYYY] =...
                GT_2DStokes (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    %singluarity tretment
    %GXX_anLin = zeros(el);
    %GYY_anLin = zeros(el);
    [SXX,SYY,QXXX,QXXY,QXYY,QYYY,GXX_anLin,GYY_anLin] =...
        Stokes2D_ST_linear(SXX,SYY,QXXX,QXXY,QXYY,QYYY,0,globalX,globalY,X0,Y0,deltaL,el);

    %INTEGRATION
    temp = reshape(repmat(deltaL/2,6,1),1,6*el);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',el,N);
    PHIA = repmat(phia',el,N);
    PHIB = repmat(phib',el,N);
            
    intSXXa = cumsum(SXX.*manyGW.*manyDelta.*PHIA);
    intSXXb = cumsum(SXX.*manyGW.*manyDelta.*PHIB);
    intSXYa = cumsum(SXY.*manyGW.*manyDelta.*PHIA);
    intSXYb = cumsum(SXY.*manyGW.*manyDelta.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*manyDelta.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*manyDelta.*PHIB);
    intQXXXa = cumsum(QXXX.*manyGW.*manyDelta.*PHIA);
    intQXXXb = cumsum(QXXX.*manyGW.*manyDelta.*PHIB);
    intQXXYa = cumsum(QXXY.*manyGW.*manyDelta.*PHIA);
    intQXXYb = cumsum(QXXY.*manyGW.*manyDelta.*PHIB);
    intQXYYa = cumsum(QXYY.*manyGW.*manyDelta.*PHIA);
    intQXYYb = cumsum(QXYY.*manyGW.*manyDelta.*PHIB);
    intQYYYa = cumsum(QYYY.*manyGW.*manyDelta.*PHIA);
    intQYYYb = cumsum(QYYY.*manyGW.*manyDelta.*PHIB);
    
    subSXXa(2:el,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el,:) = intSXYb(6:6:end-6,:);
    subSYYa(2:el,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el,:) = intQXXYb(6:6:end-6,:);
    subQXYYa(2:el,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el,:) = intQXYYb(6:6:end-6,:);
    subQYYYa(2:el,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el,:) = intQYYYb(6:6:end-6,:);
    
    %compute integral (TAKE IN ACCOUNT THAT LAST ELEMENT TOUCHES THE FIRST NODE)    
    temp = (intSXXb(6:6:end,:)-subSXXb)';
    GXX = (intSXXa(6:6:end,:)-subSXXa)' + [temp(:,end) temp(:,1:end-1)] + GXX_anLin;
    temp = (intSXYb(6:6:end,:)-subSXYb)';
    GXY = (intSXYa(6:6:end,:)-subSXYa)' + [temp(:,end) temp(:,1:end-1)];
    GYX = GXY;
    temp = (intSYYb(6:6:end,:)-subSYYb)';
    GYY = (intSYYa(6:6:end,:)-subSYYa)' + [temp(:,end) temp(:,1:end-1)] + GYY_anLin;
    
    %multiplicate normals
    R1 = repmat(r(1,:),N,1);
    R2 = repmat(r(2,:),N,1);
    
    temp = (intQXXXb(6:6:end,:)'-subQXXXb').*R1 + (intQXXYb(6:6:end,:)'-subQXXYb').*R2;
    A11 = (intQXXXa(6:6:end,:)'-subQXXXa').*R1 + (intQXXYa(6:6:end,:)'-subQXXYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    temp = (intQXXYb(6:6:end,:)'-subQXXYb').*R1 + (intQXYYb(6:6:end,:)'-subQXYYb').*R2;
    A12 = (intQXXYa(6:6:end,:)'-subQXXYa').*R1 + (intQXYYa(6:6:end,:)'-subQXYYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    temp = (intQXXYb(6:6:end,:)'-subQXXYb').*R1 + (intQXYYb(6:6:end,:)'-subQXYYb').*R2;
    A21 = (intQXXYa(6:6:end,:)'-subQXXYa').*R1 + (intQXYYa(6:6:end,:)'-subQXYYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    temp = (intQXYYb(6:6:end,:)'-subQXYYb').*R1 + (intQYYYb(6:6:end,:)'-subQYYYb').*R2;
    A22 = (intQXYYa(6:6:end,:)'-subQXYYa').*R1 + (intQYYYa(6:6:end,:)'-subQYYYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    
    %T = toc

end