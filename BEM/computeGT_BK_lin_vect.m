%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_BK_lin_vect(x,y,k)

    %tic
    
    %number of singularities
    N = numel(x)-1;
    el = N;
    
    %compute normals
    DX = x(2:end)-x(1:end-1);
    DY = y(2:end)-y(1:end-1);
    r = [DY./sqrt(DX.*DX+DY.*DY); -DX./sqrt(DX.*DX+DY.*DY)];
    
    %gauss points and weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x(1:N);
    Y0 = y(1:N);
    
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
        
    %not clear what is X0!!!
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
                GT_BK_2D (globalXmatr,globalYmatr,X0matr,Y0matr,k);
    
%     for i = 2:N-2
%        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
%        SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
%               
%        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
%        SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
%     end
%     
%     %because I don't have to treat the singularity
%     SXX(1:6,2) = SXX(1:6,2) + 2*log(sqrt((globalX(1:6)-X0(2)).^2+(globalY(1:6)-Y0(2)).^2))'-1;
%     SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))'-1;
%     
%     SYY(1:6,2) = SYY(1:6,2) + 2*log(sqrt((globalX(1:6)-X0(2)).^2+(globalY(1:6)-Y0(2)).^2))'-1;
%     SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))'-1;

    temp = reshape(repmat(deltaL/2,6,1),1,6*el);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',el,N);
    PHIA = repmat(phia',el,N);
    PHIB = repmat(phib',el,N);
            
    intSXXa = cumsum(SXX.*manyGW.*manyDelta.*PHIA);
    intSXXb = cumsum(SXX.*manyGW.*manyDelta.*PHIB);
    intSXYa = cumsum(SXY.*manyGW.*manyDelta.*PHIA);
    intSXYb = cumsum(SXY.*manyGW.*manyDelta.*PHIB);
    intSYXa = cumsum(SYX.*manyGW.*manyDelta.*PHIA);
    intSYXb = cumsum(SYX.*manyGW.*manyDelta.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*manyDelta.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*manyDelta.*PHIB);
    intQXXXa = cumsum(QXXX.*manyGW.*manyDelta.*PHIA);
    intQXXXb = cumsum(QXXX.*manyGW.*manyDelta.*PHIB);
    intQXXYa = cumsum(QXXY.*manyGW.*manyDelta.*PHIA);
    intQXXYb = cumsum(QXXY.*manyGW.*manyDelta.*PHIB);
    intQXYXa = cumsum(QXYX.*manyGW.*manyDelta.*PHIA);
    intQXYXb = cumsum(QXYX.*manyGW.*manyDelta.*PHIB);
    intQXYYa = cumsum(QXYY.*manyGW.*manyDelta.*PHIA);
    intQXYYb = cumsum(QXYY.*manyGW.*manyDelta.*PHIB);
    intQYXXa = cumsum(QYXX.*manyGW.*manyDelta.*PHIA);
    intQYXXb = cumsum(QYXX.*manyGW.*manyDelta.*PHIB);
    intQYXYa = cumsum(QYXY.*manyGW.*manyDelta.*PHIA);
    intQYXYb = cumsum(QYXY.*manyGW.*manyDelta.*PHIB);
    intQYYXa = cumsum(QYYX.*manyGW.*manyDelta.*PHIA);
    intQYYXb = cumsum(QYYX.*manyGW.*manyDelta.*PHIB);
    intQYYYa = cumsum(QYYY.*manyGW.*manyDelta.*PHIA);
    intQYYYb = cumsum(QYYY.*manyGW.*manyDelta.*PHIB);
    
    subSXXa(2:el,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:el,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:el,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:el,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el,:) = intQXXYb(6:6:end-6,:);
    subQXYXa(2:el,:) = intQXYXa(6:6:end-6,:);
    subQXYXb(2:el,:) = intQXYXb(6:6:end-6,:);
    subQXYYa(2:el,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el,:) = intQXYYb(6:6:end-6,:);
    subQYXXa(2:el,:) = intQYXXa(6:6:end-6,:);
    subQYXXb(2:el,:) = intQYXXb(6:6:end-6,:);
    subQYXYa(2:el,:) = intQYXYa(6:6:end-6,:);
    subQYXYb(2:el,:) = intQYXYb(6:6:end-6,:);
    subQYYXa(2:el,:) = intQYYXa(6:6:end-6,:);
    subQYYXb(2:el,:) = intQYYXb(6:6:end-6,:);
    subQYYYa(2:el,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el,:) = intQYYYb(6:6:end-6,:);
    
    %compute integral
    GXX = (intSXXa(6:6:end,:)-subSXXa)';
    GXX(:,[2:end 1]) = GXX(:,[2:end 1]) + intSXXb(6:6:end,:)'-subSXXb';
    GXY = intSXYa(6:6:end,:)'-subSXYa';
    GXY(:,[2:end 1]) = GXY(:,[2:end 1]) + intSXYb(6:6:end,:)'-subSXYb';
    GYX = intSYXa(6:6:end,:)'-subSYXa';
    GYX(:,[2:end 1]) = GYX(:,[2:end 1]) + intSYXb(6:6:end,:)'-subSYXb';
    GYY = intSYYa(6:6:end,:)'-subSYYa';
    GYY(:,[2:end 1]) = GYY(:,[2:end 1]) + intSYYb(6:6:end,:)'-subSYYb';
    
    %moltiplicate normals
    R1 = repmat(r(1,:),N,1);
    R2 = repmat(r(2,:),N,1);
    
    A11 = (intQXXXa(6:6:end,:)'-subQXXXa').*R1 + (intQXXYa(6:6:end,:)'-subQXXYa').*R2;
    A11([end 1:end-1]) = A11 + (intQXXXb(6:6:end,:)'-subQXXXb').*R1 + (intQXXYb(6:6:end,:)'-subQXXYb').*R2;
    A12 = (intQXYXa(6:6:end,:)'-subQXYXa').*R1 + (intQXYYa(6:6:end,:)'-subQXYYa').*R2;
    A12([end 1:end-1]) = A12 + (intQXYXb(6:6:end,:)'-subQXYXb').*R1 + (intQXYYb(6:6:end,:)'-subQXYYb').*R2;
    A21 = (intQYXXa(6:6:end,:)'-subQYXXa').*R1 + (intQYXYa(6:6:end,:)'-subQYXYa').*R2;
    A21([end 1:end-1]) = A21 + (intQYXXb(6:6:end,:)'-subQYXXb').*R1 + (intQYXYb(6:6:end,:)'-subQYXYb').*R2;
    A22 = (intQYYXa(6:6:end,:)'-subQYYXa').*R1 + (intQYYYa(6:6:end,:)'-subQYYYa').*R2;
    A22([end 1:end-1]) = A22 + (intQYYXb(6:6:end,:)'-subQYYXb').*R1 + (intQYYYb(6:6:end,:)'-subQYYYb').*R2;
    
    %singularity treatment
%     GXX(2:N-2,2:N-2) = GXX(2:N-2,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
%     GXX(2:N-2,3:N-1) = GXX(2:N-2,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
%     GXX(3:N-1,2:N-2) = GXX(3:N-1,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
%     GXX(3:N-1,3:N-1) = GXX(3:N-1,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
%     
%     GYY(2:N-2,2:N-2) = GYY(2:N-2,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
%     GYY(2:N-2,3:N-1) = GYY(2:N-2,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
%     GYY(3:N-1,2:N-2) = GYY(3:N-1,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
%     GYY(3:N-1,3:N-1) = GYY(3:N-1,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
%         
%     %because singularity on the axis don't has to be treated
%     GXX(2,1) = GXX(2,1) - deltaL(1)*log(deltaL(1))+deltaL(1);
%     GXX(2,2) = GXX(2,2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
%     GXX(N-1,N-1) = GXX(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
%     GXX(N-1,N) = GXX(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
%     
%     GYY(2,1) = GYY(2,1) - deltaL(1)*log(deltaL(1))+deltaL(1);
%     GYY(2,2) = GYY(2,2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
%     GYY(N-1,N-1) = GYY(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
%     GYY(N-1,N) = GYY(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);

    
    %T = toc

end