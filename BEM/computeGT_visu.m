%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY,PXX,PXY,PYX,PYY] = computeGT_visu(x,y,X0,Y0)

    %tic

    %number of singularities
    N = numel(X0);
    %number of elements
    el = numel(x)-1;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX = repmat(GP,1,el).*reshape((repmat(deltaX/2,6,1)),1,6*el);
    GPY = repmat(GP,1,el).*reshape((repmat(deltaY/2,6,1)),1,6*el);
    
    XX = repmat((x(1:end-1)+x(2:end))/2,6,1);
    YY = repmat((y(1:end-1)+y(2:end))/2,6,1);
    globalX = XX(:)'+GPX;
    globalY = YY(:)'+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,1,6*el)';
    Y0matr = repmat(Y0,1,6*el)';
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %not clear what is X0!!!
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,PXX,PXY,PYX,PYY] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);

    temp = reshape(repmat(deltaL/2,6,1),6*el,1);
    manyDelta = repmat(temp,1,N);
    manyGW = repmat(GW',el,N);
    
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
    
    subSXX(2:el,:) = intSXX(6:6:end-6,:);
    subSXY(2:el,:) = intSXY(6:6:end-6,:);
    subSYX(2:el,:) = intSYX(6:6:end-6,:);
    subSYY(2:el,:) = intSYY(6:6:end-6,:);
    subQXXX(2:el,:) = intQXXX(6:6:end-6,:);
    subQXXY(2:el,:) = intQXXY(6:6:end-6,:);
    subQXYX(2:el,:) = intQXYX(6:6:end-6,:);
    subQXYY(2:el,:) = intQXYY(6:6:end-6,:);
    subQYXX(2:el,:) = intQYXX(6:6:end-6,:);
    subQYXY(2:el,:) = intQYXY(6:6:end-6,:);
    subQYYX(2:el,:) = intQYYX(6:6:end-6,:);
    subQYYY(2:el,:) = intQYYY(6:6:end-6,:);
    
    GXX = intSXX(6:6:end,:)'-subSXX';
    GXY = intSXY(6:6:end,:)'-subSXY';
    GYX = intSYX(6:6:end,:)'-subSYX';
    GYY = intSYY(6:6:end,:)'-subSYY';
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