%compute Green's functions and associated stress tensor everywhere is
%needed for Brinkmann

function [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_2DStokes_interface(x,y,fixed)

    %tic

    %number of singularities
    N = numel(x)-2;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = [(x(1:fixed)+x(2:fixed+1))/2 (x(fixed+2:end-1)+x(fixed+3:end))/2];
    Y0 = [(y(1:fixed)+y(2:fixed+1))/2 (y(fixed+2:end-1)+y(fixed+3:end))/2];
    
    %moltiplicate sing points coordinates
    tempx = repmat(X0,6,1);
    XX = reshape(tempx,1,6*N);
    tempy = repmat(Y0,6,1);
    YY = reshape(tempy,1,6*N);
    
    % points where I perform gauss integration
    deltaX = [x(2:fixed+1)-x(1:fixed) x(fixed+3:end)-x(fixed+2:end-1)];
    deltaY = [y(2:fixed+1)-y(1:fixed) y(fixed+3:end)-y(fixed+2:end-1)];
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX = repmat(GP,1,N).*reshape((repmat(deltaX/2,6,1)),1,6*N);
    GPY = repmat(GP,1,N).*reshape((repmat(deltaY/2,6,1)),1,6*N);
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,6*N,1);
    Y0matr = repmat(Y0,6*N,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %not clear what is X0!!!
    %[SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
    [SXX,SXY,SYY,QXXX,QXXY,QXYY,QYYY] = ...
                GT_2DStokes(globalXmatr,globalYmatr,X0matr,Y0matr);
            
    %symmetry
    SYX = SXY;
    QXYX = QXXY;
    QYXX = QXXY;
    QYXY = QXYY;
    QYYX = QXYY;
            
    for i = 1:N
        
        %singularity treatment
        %SINGLE LAYER
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i)+log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))';
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i)+log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))';
        
        %DOUBLE LAYER
        QXXX(1+6*(i-1):6*i,i) = zeros(6,1);
        QXXY(1+6*(i-1):6*i,i) = zeros(6,1);
        QXYX(1+6*(i-1):6*i,i) = zeros(6,1);
        QXYY(1+6*(i-1):6*i,i) = zeros(6,1);
        QYXX(1+6*(i-1):6*i,i) = zeros(6,1);
        QYXY(1+6*(i-1):6*i,i) = zeros(6,1);
        QYYX(1+6*(i-1):6*i,i) = zeros(6,1);
        QYYY(1+6*(i-1):6*i,i) = zeros(6,1);
    end

    temp = reshape(repmat(deltaL/2,6,1),1,6*N);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',N,N);
    
    %INTEGRATION
    intSXX = cumsum(SXX.*manyGW.*manyDelta);
    intSXY = cumsum(SXY.*manyGW.*manyDelta);
    intSYX = cumsum(SYX.*manyGW.*manyDelta);
    intSYY = cumsum(SYY.*manyGW.*manyDelta);
    intQXXX = cumsum(QXXX.*manyGW.*manyDelta);
    intQXXY = cumsum(QXXY.*manyGW.*manyDelta);
    intQXYX = cumsum(QXYX.*manyGW.*manyDelta);
    intQXYY = cumsum(QXYY.*manyGW.*manyDelta);
    intQYXX = cumsum(QYXX.*manyGW.*manyDelta);
    intQYXY = cumsum(QYXY.*manyGW.*manyDelta);
    intQYYX = cumsum(QYYX.*manyGW.*manyDelta);
    intQYYY = cumsum(QYYY.*manyGW.*manyDelta);
    
    subSXX(2:N,:) = intSXX(6:6:end-6,:);
    subSXY(2:N,:) = intSXY(6:6:end-6,:);
    subSYX(2:N,:) = intSYX(6:6:end-6,:);
    subSYY(2:N,:) = intSYY(6:6:end-6,:);
    subQXXX(2:N,:) = intQXXX(6:6:end-6,:);
    subQXXY(2:N,:) = intQXXY(6:6:end-6,:);
    subQXYX(2:N,:) = intQXYX(6:6:end-6,:);
    subQXYY(2:N,:) = intQXYY(6:6:end-6,:);
    subQYXX(2:N,:) = intQYXX(6:6:end-6,:);
    subQYXY(2:N,:) = intQYXY(6:6:end-6,:);
    subQYYX(2:N,:) = intQYYX(6:6:end-6,:);
    subQYYY(2:N,:) = intQYYY(6:6:end-6,:);
    
    GXX = intSXX(6:6:end,:)'-subSXX'+diag(-deltaL.*log(deltaL/2)+deltaL);
    GXY = intSXY(6:6:end,:)'-subSXY';
    GYX = intSYX(6:6:end,:)'-subSYX';
    GYY = intSYY(6:6:end,:)'-subSYY'+diag(-deltaL.*log(deltaL/2)+deltaL);
    TXXX = intQXXX(6:6:end,:)'-subQXXX';
    TXXY = intQXXY(6:6:end,:)'-subQXXY';
    TXYX = intQXYX(6:6:end,:)'-subQXYX';
    TXYY = intQXYY(6:6:end,:)'-subQXYY';
    TYXX = intQYXX(6:6:end,:)'-subQYXX';
    TYXY = intQYXY(6:6:end,:)'-subQYXY';
    TYYX = intQYYX(6:6:end,:)'-subQYYX';
    TYYY = intQYYY(6:6:end,:)'-subQYYY';
    
    %T = toc

end