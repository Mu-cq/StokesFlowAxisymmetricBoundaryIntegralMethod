%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_2DStokes_interface_linear(x,y,walls,r,PARAM)

    %tic

    %number of singularities and element
    N = numel(x)-2;
    el = N;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %hat functions
    phia = 1-(GP+1)/2;
    phib = (GP+1)/2;
    
    % point where the variable will be stored
    X0 = [(x(1:walls)+x(2:walls+1))/2 x(walls+2:end-1)];
    Y0 = [(y(1:walls)+y(2:walls+1))/2 y(walls+2:end-1)];
    
    % points where I perform gauss integration
    deltaX = [x(2:walls+1)-x(1:walls) x(walls+3:end)-x(walls+2:end-1)];
    deltaY = [y(2:walls+1)-y(1:walls) y(walls+3:end)-y(walls+2:end-1)];
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %moltiplicate X0 and Y0
    tempx = repmat([X0(1:walls) (x(walls+2:end-1)+x(walls+3:end))/2],6,1);
    XX = tempx(:)';
    tempy = repmat([Y0(1:walls) (y(walls+2:end-1)+y(walls+3:end))/2],6,1);
    YY = tempy(:)';
    
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
    
    if PARAM.cfunction == 0
        
    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %not clear what is X0!!!
    [SXX,SXY,SYY,QXXX,QXXY,QXYY,QYYY] =...
                GT_2DStokes (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    elseif PARAM.cfunction == 1
            
    [SXX,SXY,SYY,QXXX,QXXY,QXYY,QYYY]...
        = GT_2DStokes_cfunction(globalX,globalY,X0,Y0);
    
    SXX = reshape(SXX',numel(X0),numel(globalX))';
    SXY = reshape(SXY',numel(X0),numel(globalX))';
    SYY = reshape(SYY',numel(X0),numel(globalX))';
    QXXX = reshape(QXXX',numel(X0),numel(globalX))';
    QXXY = reshape(QXXY',numel(X0),numel(globalX))';
    QXYY = reshape(QXYY',numel(X0),numel(globalX))';
    QYYY = reshape(QYYY',numel(X0),numel(globalX))';
    
    end
    
    %singular treatment constant elements WALL
    %GXX_an = 0;
    %GYY_an = 0;
    [SXX,SYY,QXXX,QXXY,QXYY,QYYY,GXX_an,GYY_an] =...
       Stokes2D_ST_const(SXX,SYY,QXXX,QXXY,QXYY,QYYY,walls,globalX,globalY,X0,Y0,deltaL,el);

    %singular treatment constant elements WALL
    %GXX_anLin = zeros(N-walls);
    %GYY_anLin = zeros(N-walls);
    [SXX,SYY,QXXX,QXXY,QXYY,QYYY,GXX_anLin,GYY_anLin] = ...
        Stokes2D_ST_linear(SXX,SYY,QXXX,QXXY,QXYY,QYYY,walls,globalX,globalY,X0,Y0,deltaL(walls+1:end),el);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT ELEMENTS INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    temp = reshape(repmat(deltaL(1:walls)/2,6,1),1,6*walls);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',walls,N);
    
    intSXX = cumsum(SXX(1:6*walls,:).*manyGW.*manyDelta);
    intSXY = cumsum(SXY(1:6*walls,:).*manyGW.*manyDelta);
    intSYY = cumsum(SYY(1:6*walls,:).*manyGW.*manyDelta);
    intQXXX = cumsum(QXXX(1:6*walls,:).*manyGW.*manyDelta);
    intQXXY = cumsum(QXXY(1:6*walls,:).*manyGW.*manyDelta);
    intQXYY = cumsum(QXYY(1:6*walls,:).*manyGW.*manyDelta);
    intQYYY = cumsum(QYYY(1:6*walls,:).*manyGW.*manyDelta);
    
    subSXX(2:walls,:) = intSXX(6:6:end-6,:);
    subSXY(2:walls,:) = intSXY(6:6:end-6,:);
    subSYY(2:walls,:) = intSYY(6:6:end-6,:);
    subQXXX(2:walls,:) = intQXXX(6:6:end-6,:);
    subQXXY(2:walls,:) = intQXXY(6:6:end-6,:);
    subQXYY(2:walls,:) = intQXYY(6:6:end-6,:);
    subQYYY(2:walls,:) = intQYYY(6:6:end-6,:);
    
    GXX(:,1:walls) = intSXX(6:6:end,:)'-subSXX' + GXX_an;
    GXY(:,1:walls) = intSXY(6:6:end,:)'-subSXY';
    GYY(:,1:walls) = intSYY(6:6:end,:)'-subSYY' + GYY_an;
    TXXX(:,1:walls) = intQXXX(6:6:end,:)'-subQXXX';
    TXXY(:,1:walls) = intQXXY(6:6:end,:)'-subQXXY';
    TXYY(:,1:walls) = intQXYY(6:6:end,:)'-subQXYY';
    TYYY(:,1:walls) = intQYYY(6:6:end,:)'-subQYYY';
    
    R1 = repmat(r(1,1:walls),N,1);
    R2 = repmat(r(2,1:walls),N,1);
    
    A11(:,1:walls) = TXXX.*R1+TXXY.*R2;
    A12(:,1:walls) = TXXY.*R1+TXYY.*R2;
    A21(:,1:walls) = TXXY.*R1+TXYY.*R2;
    A22(:,1:walls) = TXYY.*R1+TYYY.*R2;
    
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
    intSYYa = cumsum(SYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intSYYb = cumsum(SYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXXXa = cumsum(QXXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXXXb = cumsum(QXXX(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXXYa = cumsum(QXXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXXYb = cumsum(QXXY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQXYYa = cumsum(QXYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQXYYb = cumsum(QXYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    intQYYYa = cumsum(QYYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIA);
    intQYYYb = cumsum(QYYY(1+6*walls:end,:).*manyGW.*manyDelta.*PHIB);
    
    subSXXa(2:el-walls,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el-walls,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el-walls,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el-walls,:) = intSXYb(6:6:end-6,:);
    subSYYa(2:el-walls,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el-walls,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el-walls,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el-walls,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el-walls,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el-walls,:) = intQXXYb(6:6:end-6,:);
    subQXYYa(2:el-walls,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el-walls,:) = intQXYYb(6:6:end-6,:);
    subQYYYa(2:el-walls,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el-walls,:) = intQYYYb(6:6:end-6,:);
        
    %compute integral (TAKE IN ACCOUNT THAT LAST ELEMENT TOUCHES THE FIRST NODE)    
    temp = (intSXXb(6:6:end,:)-subSXXb)';
    GXX(:,walls+1:end) = (intSXXa(6:6:end,:)-subSXXa)' + [temp(:,end) temp(:,1:end-1)] + ...
        [zeros(walls,N-walls); GXX_anLin];
    temp = (intSXYb(6:6:end,:)-subSXYb)';
    GXY(:,walls+1:end) = (intSXYa(6:6:end,:)-subSXYa)' + [temp(:,end) temp(:,1:end-1)];
    GYX = GXY;
    temp = (intSYYb(6:6:end,:)-subSYYb)';
    GYY(:,walls+1:end) = (intSYYa(6:6:end,:)-subSYYa)' + [temp(:,end) temp(:,1:end-1)] + ...
        [zeros(walls,N-walls); GYY_anLin];
    
    %multiplicate normals
    R1 = repmat(r(1,walls+1:end),N,1);
    R2 = repmat(r(2,walls+1:end),N,1);
    
    temp = (intQXXXb(6:6:end,:)'-subQXXXb').*R1 + (intQXXYb(6:6:end,:)'-subQXXYb').*R2;
    A11(:,walls+1:end) = (intQXXXa(6:6:end,:)'-subQXXXa').*R1 + (intQXXYa(6:6:end,:)'-subQXXYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    temp = (intQXXYb(6:6:end,:)'-subQXXYb').*R1 + (intQXYYb(6:6:end,:)'-subQXYYb').*R2;
    A12(:,walls+1:end) = (intQXXYa(6:6:end,:)'-subQXXYa').*R1 + (intQXYYa(6:6:end,:)'-subQXYYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    temp = (intQXXYb(6:6:end,:)'-subQXXYb').*R1 + (intQXYYb(6:6:end,:)'-subQXYYb').*R2;
    A21(:,walls+1:end) = (intQXXYa(6:6:end,:)'-subQXXYa').*R1 + (intQXYYa(6:6:end,:)'-subQXYYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    temp = (intQXYYb(6:6:end,:)'-subQXYYb').*R1 + (intQYYYb(6:6:end,:)'-subQYYYb').*R2;
    A22(:,walls+1:end) = (intQXYYa(6:6:end,:)'-subQXYYa').*R1 + (intQYYYa(6:6:end,:)'-subQYYYa').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %T = toc

end