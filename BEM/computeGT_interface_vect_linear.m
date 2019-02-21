%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2,Iaxis] = computeGT_interface_vect_linear(x,y,walls,r)

    %tic

    %number of singularities and element
    N = numel(x)-1;
    el = N-1;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    T1 = zeros(N,el);
    T2 = zeros(N,el);
    D1 = zeros(N,el);
    D2 = zeros(N,el);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %hat functions
    phia = 1-(GP+1)/2;
    phib = (GP+1)/2;
    
    % point where the variable will be stored
    X0 = [(x(1:walls)+x(2:walls+1))/2 x(walls+2:end)];
    Y0 = [(y(1:walls)+y(2:walls+1))/2 y(walls+2:end)];
    
    %moltiplicate X0 and Y0
    tempx = repmat([X0(1:walls) (X0(walls+1:end-1)+X0(walls+2:end))/2],6,1);
    XX = reshape(tempx,1,6*el);
    tempy = repmat([Y0(1:walls) (Y0(walls+1:end-1)+Y0(walls+2:end))/2],6,1);
    YY = reshape(tempy,1,6*el);
    
    % points where I perform gauss integration
    deltaX = [x(2:walls+1)-x(1:walls) x(walls+3:end)-x(walls+2:end-1)];
    deltaY = [y(2:walls+1)-y(1:walls) y(walls+3:end)-y(walls+2:end-1)];
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX = repmat(GP,1,el).*reshape((repmat(deltaX/2,6,1)),1,6*el);
    GPY = repmat(GP,1,el).*reshape((repmat(deltaY/2,6,1)),1,6*el);
    
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
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,PXX,PXY,PYX,PYY ,Iaxis] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
        
    % singularity treatment const elem
    for k = 1:walls
        
        SXX(1+6*(k-1):6*k,k) = SXX(1+6*(k-1):6*k,k)...
            +2*log(sqrt((globalX(1+6*(k-1):6*k)-X0(k)).^2+(globalY(1+6*(k-1):6*k)-Y0(k)).^2)')-1;
        SYY(1+6*(k-1):6*k,k) = SYY(1+6*(k-1):6*k,k)...
            +2*log(sqrt((globalX(1+6*(k-1):6*k)-X0(k)).^2+(globalY(1+6*(k-1):6*k)-Y0(k)).^2)')-1;
        
    end
    
    % singularity treatment lin elem
    for i = walls+2:el-1
        
       SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) ...
            + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) ...
            + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
              
       SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) ...
            + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) ...
            + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
        
    end
     
    %because I don't have to treat the singularity ON THE AXIS
    SXX(6*walls+1:6*walls+6,walls+2) = SXX(6*walls+1:6*walls+6,walls+2) ...
        + 2*log(sqrt((globalX(6*walls+1:6*walls+6)-X0(walls+2)).^2+(globalY(6*walls+1:6*walls+6)-Y0(walls+2)).^2))'-1;
    SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) ...
        + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))'-1;
    
    SYY(6*walls+1:6*walls+6,walls+2) = SYY(6*walls+1:6*walls+6,walls+2) ...
        + 2*log(sqrt((globalX(6*walls+1:6*walls+6)-X0(walls+2)).^2+(globalY(6*walls+1:6*walls+6)-Y0(walls+2)).^2))'-1;
    SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) ...
        + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))'-1;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT ELEMENTS INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    temp = reshape(repmat(deltaL(1:walls)/2,6,1),1,6*walls);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',walls,N);
    
    intSXX = cumsum(SXX(1:6*walls,:).*manyGW.*manyDelta);
    intSXY = cumsum(SXY(1:6*walls,:).*manyGW.*manyDelta);
    intSYX = cumsum(SYX(1:6*walls,:).*manyGW.*manyDelta);
    intSYY = cumsum(SYY(1:6*walls,:).*manyGW.*manyDelta);
    intQXXX = cumsum(QXXX(1:6*walls,:).*manyGW.*manyDelta);
    intQXXY = cumsum(QXXY(1:6*walls,:).*manyGW.*manyDelta);
    intQXYX = cumsum(QXYX(1:6*walls,:).*manyGW.*manyDelta);
    intQXYY = cumsum(QXYY(1:6*walls,:).*manyGW.*manyDelta);
    intQYXX = cumsum(QYXX(1:6*walls,:).*manyGW.*manyDelta);
    intQYXY = cumsum(QYXY(1:6*walls,:).*manyGW.*manyDelta);
    intQYYX = cumsum(QYYX(1:6*walls,:).*manyGW.*manyDelta);
    intQYYY = cumsum(QYYY(1:6*walls,:).*manyGW.*manyDelta);
    intPXX = cumsum(PXX(1:6*walls,:).*manyGW.*manyDelta);
    intPXY = cumsum(PXY(1:6*walls,:).*manyGW.*manyDelta);
    intPYX = cumsum(PYX(1:6*walls,:).*manyGW.*manyDelta);
    intPYY = cumsum(PYY(1:6*walls,:).*manyGW.*manyDelta);
    
    subSXX(2:walls,:) = intSXX(6:6:end-6,:);
    subSXY(2:walls,:) = intSXY(6:6:end-6,:);
    subSYX(2:walls,:) = intSYX(6:6:end-6,:);
    subSYY(2:walls,:) = intSYY(6:6:end-6,:);
    subQXXX(2:walls,:) = intQXXX(6:6:end-6,:);
    subQXXY(2:walls,:) = intQXXY(6:6:end-6,:);
    subQXYX(2:walls,:) = intQXYX(6:6:end-6,:);
    subQXYY(2:walls,:) = intQXYY(6:6:end-6,:);
    subQYXX(2:walls,:) = intQYXX(6:6:end-6,:);
    subQYXY(2:walls,:) = intQYXY(6:6:end-6,:);
    subQYYX(2:walls,:) = intQYYX(6:6:end-6,:);
    subQYYY(2:walls,:) = intQYYY(6:6:end-6,:);
    subPXX(2:walls,:) = intPXX(6:6:end-6,:);
    subPXY(2:walls,:) = intPXY(6:6:end-6,:);
    subPYX(2:walls,:) = intPYX(6:6:end-6,:);
    subPYY(2:walls,:) = intPYY(6:6:end-6,:);
    
    GXX(:,1:walls) = intSXX(6:6:end,:)'-subSXX' ...
        + [diag(2*(-2*deltaL(1:walls)/2.*log(deltaL(1:walls)/2)+3*deltaL(1:walls)/2)); zeros(el-walls+1,walls)];
    GXY(:,1:walls) = intSXY(6:6:end,:)'-subSXY';
    GYX(:,1:walls) = intSYX(6:6:end,:)'-subSYX';
    GYY(:,1:walls) = intSYY(6:6:end,:)'-subSYY' ...
        + [diag(2*(-2*deltaL(1:walls)/2.*log(deltaL(1:walls)/2)+3*deltaL(1:walls)/2)); zeros(el-walls+1,walls)];
    TXXX(:,1:walls) = intQXXX(6:6:end,:)'-subQXXX';
    TXXY(:,1:walls) = intQXXY(6:6:end,:)'-subQXXY';
    TXYX(:,1:walls) = intQXYX(6:6:end,:)'-subQXYX';
    TXYY(:,1:walls) = intQXYY(6:6:end,:)'-subQXYY';
    TYXX(:,1:walls) = intQYXX(6:6:end,:)'-subQYXX';
    TYXY(:,1:walls) = intQYXY(6:6:end,:)'-subQYXY';
    TYYX(:,1:walls) = intQYYX(6:6:end,:)'-subQYYX';
    TYYY(:,1:walls) = intQYYY(6:6:end,:)'-subQYYY';
    DXX(:,1:walls) = intPXX(6:6:end,:)'-subPXX';
    DXY(:,1:walls) = intPXY(6:6:end,:)'-subPXY';
    DYX(:,1:walls) = intPYX(6:6:end,:)'-subPYX';
    DYY(:,1:walls) = intPYY(6:6:end,:)'-subPYY';
    
    R1 = repmat(r(1,1:walls),N,1);
    R2 = repmat(r(2,1:walls),N,1);
    
    A11(:,1:walls) = TXXX.*R1+TXXY.*R2;
    A12(:,1:walls) = TXYX.*R1+TXYY.*R2;
    A21(:,1:walls) = TYXX.*R1+TYXY.*R2;
    A22(:,1:walls) = TYYX.*R1+TYYY.*R2;
    T1(:,1:walls) = A11(:,1:walls);
    T2(:,1:walls) = A21(:,1:walls);
    D1(:,1:walls) = DXX.*R1+DXY.*R2;
    D2(:,1:walls) = DYX.*R1+DYY.*R2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR ELEMENTS INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    temp = reshape(repmat(deltaL(walls+1:end)/2,6,1),1,6*(el-walls));
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',el-walls,N);
    PHIA = repmat(phia',el-walls,N);
    PHIB = repmat(phib',el-walls,N);
            
    intSXXa = cumsum(SXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intSXXb = cumsum(SXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intSXYa = cumsum(SXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intSXYb = cumsum(SXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intSYXa = cumsum(SYX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intSYXb = cumsum(SYX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intSYYa = cumsum(SYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intSYYb = cumsum(SYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXXXa = cumsum(QXXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXXXb = cumsum(QXXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXXYa = cumsum(QXXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXXYb = cumsum(QXXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXYXa = cumsum(QXYX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXYXb = cumsum(QXYX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXYYa = cumsum(QXYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXYYb = cumsum(QXYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQYXXa = cumsum(QYXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQYXXb = cumsum(QYXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQYXYa = cumsum(QYXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQYXYb = cumsum(QYXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQYYXa = cumsum(QYYX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQYYXb = cumsum(QYYX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQYYYa = cumsum(QYYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQYYYb = cumsum(QYYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXXX_2nd = cumsum(QXXX(1+6*walls:end,:).*manyGW.*manyDelta);
    intQXXY_2nd = cumsum(QXXY(1+6*walls:end,:).*manyGW.*manyDelta);
    intQYXX_2nd = cumsum(QYXX(1+6*walls:end,:).*manyGW.*manyDelta);
    intQYXY_2nd = cumsum(QYXY(1+6*walls:end,:).*manyGW.*manyDelta);
    intPXX_2nd = cumsum(PXX(1+6*walls:end,:).*manyGW.*manyDelta);
    intPXY_2nd = cumsum(PXY(1+6*walls:end,:).*manyGW.*manyDelta);
    intPYX_2nd = cumsum(PYX(1+6*walls:end,:).*manyGW.*manyDelta);
    intPYY_2nd = cumsum(PYY(1+6*walls:end,:).*manyGW.*manyDelta);
    
    subSXXa(2:el-walls,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el-walls,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el-walls,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el-walls,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:el-walls,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:el-walls,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:el-walls,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el-walls,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el-walls,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el-walls,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el-walls,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el-walls,:) = intQXXYb(6:6:end-6,:);
    subQXYXa(2:el-walls,:) = intQXYXa(6:6:end-6,:);
    subQXYXb(2:el-walls,:) = intQXYXb(6:6:end-6,:);
    subQXYYa(2:el-walls,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el-walls,:) = intQXYYb(6:6:end-6,:);
    subQYXXa(2:el-walls,:) = intQYXXa(6:6:end-6,:);
    subQYXXb(2:el-walls,:) = intQYXXb(6:6:end-6,:);
    subQYXYa(2:el-walls,:) = intQYXYa(6:6:end-6,:);
    subQYXYb(2:el-walls,:) = intQYXYb(6:6:end-6,:);
    subQYYXa(2:el-walls,:) = intQYYXa(6:6:end-6,:);
    subQYYXb(2:el-walls,:) = intQYYXb(6:6:end-6,:);
    subQYYYa(2:el-walls,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el-walls,:) = intQYYYb(6:6:end-6,:);
    subQXXX_2nd(2:el-walls,:) = intQXXX_2nd(6:6:end-6,:);
    subQXXY_2nd(2:el-walls,:) = intQXXY_2nd(6:6:end-6,:);
    subQYXX_2nd(2:el-walls,:) = intQYXX_2nd(6:6:end-6,:);
    subQYXY_2nd(2:el-walls,:) = intQYXY_2nd(6:6:end-6,:);
    subPXX_2nd(2:el-walls,:) = intPXX_2nd(6:6:end-6,:);
    subPXY_2nd(2:el-walls,:) = intPXY_2nd(6:6:end-6,:);
    subPYX_2nd(2:el-walls,:) = intPYX_2nd(6:6:end-6,:);
    subPYY_2nd(2:el-walls,:) = intPYY_2nd(6:6:end-6,:);
    
    %compute integral
    GXX(:,walls+1:N-1) = (intSXXa(6:6:end,:)-subSXXa)';
    GXX(:,walls+2:N) = GXX(:,walls+2:N) + intSXXb(6:6:end,:)'-subSXXb';
    GXY(:,walls+1:N-1) = intSXYa(6:6:end,:)'-subSXYa';
    GXY(:,walls+2:N) = GXY(:,walls+2:N) + intSXYb(6:6:end,:)'-subSXYb';
    GYX(:,walls+1:N-1) = intSYXa(6:6:end,:)'-subSYXa';
    GYX(:,walls+2:N) = GYX(:,walls+2:N) + intSYXb(6:6:end,:)'-subSYXb';
    GYY(:,walls+1:N-1) = intSYYa(6:6:end,:)'-subSYYa';
    GYY(:,walls+2:N) = GYY(:,walls+2:N) + intSYYb(6:6:end,:)'-subSYYb';
    
    %multiplicate normals
    R1 = repmat(r(1,walls+1:end),N,1);
    R2 = repmat(r(2,walls+1:end),N,1);
    
    A11(:,walls+1:N-1) = (intQXXXa(6:6:end,:)'-subQXXXa').*R1 + (intQXXYa(6:6:end,:)'-subQXXYa').*R2;
    A11(:,walls+2:N) = A11(:,walls+2:N) + (intQXXXb(6:6:end,:)'-subQXXXb').*R1 + (intQXXYb(6:6:end,:)'-subQXXYb').*R2;
    A12(:,walls+1:N-1) = (intQXYXa(6:6:end,:)'-subQXYXa').*R1 + (intQXYYa(6:6:end,:)'-subQXYYa').*R2;
    A12(:,walls+2:N) = A12(:,walls+2:N) + (intQXYXb(6:6:end,:)'-subQXYXb').*R1 + (intQXYYb(6:6:end,:)'-subQXYYb').*R2;
    A21(:,walls+1:N-1) = (intQYXXa(6:6:end,:)'-subQYXXa').*R1 + (intQYXYa(6:6:end,:)'-subQYXYa').*R2;
    A21(:,walls+2:N) = A21(:,walls+2:N) + (intQYXXb(6:6:end,:)'-subQYXXb').*R1 + (intQYXYb(6:6:end,:)'-subQYXYb').*R2;
    A22(:,walls+1:N-1) = (intQYYXa(6:6:end,:)'-subQYYXa').*R1 + (intQYYYa(6:6:end,:)'-subQYYYa').*R2;
    A22(:,walls+2:N) = A22(:,walls+2:N) + (intQYYXb(6:6:end,:)'-subQYYXb').*R1 + (intQYYYb(6:6:end,:)'-subQYYYb').*R2;
    
    T1(:,walls+1:end) = (intQXXX_2nd(6:6:end,:)'-subQXXX_2nd').*R1 + (intQXXY_2nd(6:6:end,:)'-subQXXY_2nd').*R2;
    T2(:,walls+1:end) = (intQYXX_2nd(6:6:end,:)'-subQYXX_2nd').*R1 + (intQYXY_2nd(6:6:end,:)'-subQYXY_2nd').*R2;
    
    D1(:,walls+1:end) = (intPXX_2nd(6:6:end,:)'-subPXX_2nd').*R1 + (intPXY_2nd(6:6:end,:)'-subPXY_2nd').*R2;
    D2(:,walls+1:end) = (intPYX_2nd(6:6:end,:)'-subPYX_2nd').*R1 + (intPYY_2nd(6:6:end,:)'-subPYY_2nd').*R2;
    
    %singularity treatment
    GXX(walls+2:N-2,walls+2:N-2) = GXX(walls+2:N-2,walls+2:N-2) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+2*deltaL(walls+2:N-2));
    GXX(walls+2:N-2,walls+3:N-1) = GXX(walls+2:N-2,walls+3:N-1) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+deltaL(walls+2:N-2));
    GXX(walls+3:N-1,walls+2:N-2) = GXX(walls+3:N-1,walls+2:N-2) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+deltaL(walls+2:N-2));
    GXX(walls+3:N-1,walls+3:N-1) = GXX(walls+3:N-1,walls+3:N-1) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+2*deltaL(walls+2:N-2));
    
    GYY(walls+2:N-2,walls+2:N-2) = GYY(walls+2:N-2,walls+2:N-2) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+2*deltaL(walls+2:N-2));
    GYY(walls+2:N-2,walls+3:N-1) = GYY(walls+2:N-2,walls+3:N-1) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+deltaL(walls+2:N-2));
    GYY(walls+3:N-1,walls+2:N-2) = GYY(walls+3:N-1,walls+2:N-2)...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+deltaL(walls+2:N-2));
    GYY(walls+3:N-1,walls+3:N-1) = GYY(walls+3:N-1,walls+3:N-1) ...
        + diag(-deltaL(walls+2:N-2).*log(deltaL(walls+2:N-2))+2*deltaL(walls+2:N-2));
    
    %because singularity on the axis don't has to be treated
    GXX(walls+2,walls+1) = GXX(walls+2,walls+1) - deltaL(walls+1)*log(deltaL(walls+1))+deltaL(walls+1);
    GXX(walls+2,walls+2) = GXX(walls+2,walls+2) - deltaL(walls+1)*log(deltaL(walls+1))+2*deltaL(walls+1);
    GXX(N-1,N-1) = GXX(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
    GXX(N-1,N) = GXX(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
    
    GYY(walls+2,walls+1) = GYY(walls+2,walls+1) - deltaL(walls+1)*log(deltaL(walls+1))+deltaL(walls+1);
    GYY(walls+2,walls+2) = GYY(walls+2,walls+2) - deltaL(walls+1)*log(deltaL(walls+1))+2*deltaL(walls+1);
    GYY(N-1,N-1) = GYY(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
    GYY(N-1,N) = GYY(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %T = toc

end