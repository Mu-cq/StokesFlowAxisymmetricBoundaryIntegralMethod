%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_lin_vect2(x,y)

    %tic
    
    %number of singularities
    N = numel(x);
    el = N-1;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    
    %compute normals
    DX = x(2:end)-x(1:end-1);
    DY = y(2:end)-y(1:end-1);
    
    %normal vector
    r = [DY./sqrt(DX.*DX+DY.*DY); -DX./sqrt(DX.*DX+DY.*DY)];
    
    %multiplicate normals
    nx = r(1,:);    nx = repmat(nx,6,1);    nx = nx(:);
    ny = r(2,:);    ny = repmat(ny,6,1);    ny = ny(:);
    R1 = repmat(nx,1,N);
    R2 = repmat(ny,1,N);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x;
    Y0 = y;
    
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
    
    globalX = XX(1:end-6)+GPX;
    globalY = YY(1:end-6)+GPY;
    
    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %not clear what is X0!!!
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %double layer kernel
    Q11 = QXXX.*R1 + QXXY.*R2;
    Q12 = QXYX.*R1 + QXYY.*R2;
    Q21 = QYXX.*R1 + QYXY.*R2;
    Q22 = QYYX.*R1 + QYYY.*R2;
            
    %singularity treatment
    for i = 2:N-2
       SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
              
       SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
    end
    
    %because I don't have to treat the singularity ON THE AXIS
    SXX(1:6,2) = SXX(1:6,2) + 2*log(sqrt((globalX(1:6)-X0(2)).^2+(globalY(1:6)-Y0(2)).^2))'-1;
    SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) ...
        + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))'-1;
    
    SYY(1:6,2) = SYY(1:6,2) + 2*log(sqrt((globalX(1:6)-X0(2)).^2+(globalY(1:6)-Y0(2)).^2))'-1;
    SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) ...
        + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))'-1;

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
    intSYXa = cumsum(SYX.*manyGW.*manyDelta.*PHIA);
    intSYXb = cumsum(SYX.*manyGW.*manyDelta.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*manyDelta.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*manyDelta.*PHIB);
    intQ11a = cumsum(Q11.*manyGW.*manyDelta.*PHIA);
    intQ11b = cumsum(Q11.*manyGW.*manyDelta.*PHIB);
    intQ12a = cumsum(Q12.*manyGW.*manyDelta.*PHIA);
    intQ12b = cumsum(Q12.*manyGW.*manyDelta.*PHIB);
    intQ21a = cumsum(Q21.*manyGW.*manyDelta.*PHIA);
    intQ21b = cumsum(Q21.*manyGW.*manyDelta.*PHIB);
    intQ22a = cumsum(Q22.*manyGW.*manyDelta.*PHIA);
    intQ22b = cumsum(Q22.*manyGW.*manyDelta.*PHIB);
    
    subSXXa(2:el,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:el,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:el,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:el,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el,:) = intSYYb(6:6:end-6,:);
    subQ11a(2:el,:) = intQ11a(6:6:end-6,:);
    subQ11b(2:el,:) = intQ11b(6:6:end-6,:);
    subQ12a(2:el,:) = intQ12a(6:6:end-6,:);
    subQ12b(2:el,:) = intQ12b(6:6:end-6,:);
    subQ21a(2:el,:) = intQ21a(6:6:end-6,:);
    subQ21b(2:el,:) = intQ21b(6:6:end-6,:);
    subQ22a(2:el,:) = intQ22a(6:6:end-6,:);
    subQ22b(2:el,:) = intQ22b(6:6:end-6,:);
    
    %compute integral
    GXX(1:N,1:N-1) = (intSXXa(6:6:end,:)-subSXXa)'; %+ vertcat(zeros(1,N-1),diag(-deltaL(1:N-1).*log(deltaL(1:N-1))+deltaL(1:N-1)));
    GXX(1:N,2:N) = GXX(1:N,2:N) + intSXXb(6:6:end,:)'-subSXXb';
    GXY(1:N,1:N-1) = intSXYa(6:6:end,:)'-subSXYa';
    GXY(1:N,2:N) = GXY(1:N,2:N) + intSXYb(6:6:end,:)'-subSXYb';
    GYX(1:N,1:N-1) = intSYXa(6:6:end,:)'-subSYXa';
    GYX(1:N,2:N) = GYX(1:N,2:N) + intSYXb(6:6:end,:)'-subSYXb';
    GYY(1:N,1:N-1) = intSYYa(6:6:end,:)'-subSYYa';
    GYY(1:N,2:N) = GYY(1:N,2:N) + intSYYb(6:6:end,:)'-subSYYb';
    
    A11(:,1:N-1) = (intQ11a(6:6:end,:)-subQ11a)';
    A11(:,2:N) = A11(:,2:N) + (intQ11b(6:6:end,:)-subQ11b)';
    A12(:,1:N-1) = (intQ12a(6:6:end,:)-subQ12a)';
    A12(:,2:N) = A12(:,2:N) + (intQ12b(6:6:end,:)-subQ12b)';
    A21(:,1:N-1) = (intQ21a(6:6:end,:)-subQ21a)';
    A21(:,2:N) = A21(:,2:N) + (intQ21b(6:6:end,:)-subQ21b)';
    A22(:,1:N-1) = (intQ22a(6:6:end,:)-subQ22a)';
    A22(:,2:N) = A22(:,2:N) + (intQ22b(6:6:end,:)-subQ22b)';
    
    %singularity treatment
    GXX(2:N-2,2:N-2) = GXX(2:N-2,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
    GXX(2:N-2,3:N-1) = GXX(2:N-2,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
    GXX(3:N-1,2:N-2) = GXX(3:N-1,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
    GXX(3:N-1,3:N-1) = GXX(3:N-1,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
    
    GYY(2:N-2,2:N-2) = GYY(2:N-2,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
    GYY(2:N-2,3:N-1) = GYY(2:N-2,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
    GYY(3:N-1,2:N-2) = GYY(3:N-1,2:N-2) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+deltaL(2:N-2));
    GYY(3:N-1,3:N-1) = GYY(3:N-1,3:N-1) + diag(-deltaL(2:N-2).*log(deltaL(2:N-2))+2*deltaL(2:N-2));
        
    %because singularity on the axis don't has to be treated
    GXX(2,1) = GXX(2,1) - deltaL(1)*log(deltaL(1))+deltaL(1);
    GXX(2,2) = GXX(2,2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
    GXX(N-1,N-1) = GXX(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
    GXX(N-1,N) = GXX(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
    
    GYY(2,1) = GYY(2,1) - deltaL(1)*log(deltaL(1))+deltaL(1);
    GYY(2,2) = GYY(2,2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
    GYY(N-1,N-1) = GYY(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
    GYY(N-1,N) = GYY(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);

    
    %T = toc

end