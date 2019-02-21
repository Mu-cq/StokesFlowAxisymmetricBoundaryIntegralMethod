%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,D1,D2,Iaxis] = computeGT_interface_vect(x,y,walls,r)

    %tic

    %number of singularities
    N = numel(x)-2;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %GW = GW/2;
    
    % point where the variable will be stored
    X0 = [(x(1:walls)+x(2:walls+1))/2 (x(walls+2:end-1)+x(walls+3:end))/2];
    Y0 = [(y(1:walls)+y(2:walls+1))/2 (y(walls+2:end-1)+y(walls+3:end))/2];
    
    %moltiplicate X0 and Y0
    tempx = repmat(X0,6,1);
    XX = reshape(tempx,1,6*N);
    tempy = repmat(Y0,6,1);
    YY = reshape(tempy,1,6*N);
    
    % points where I perform gauss integration
    deltaX = [x(2:walls+1)-x(1:walls) x(walls+3:end)-x(walls+2:end-1)];
    deltaY = [y(2:walls+1)-y(1:walls) y(walls+3:end)-y(walls+2:end-1)];
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
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,PXX,PXY,PYX,PYY ,Iaxis] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
        
    for i = 1:N        
        %singularity treatment
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
    end
    
%     tempx = reshape(globalX,N,6);
%     tempy = reshape(globalY,N,6);
% 
%     SXX = SXX + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
%     SYY = SYY + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
    
    temp = reshape(repmat(deltaL/2,6,1),1,6*N);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',N,N);
    
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
    intPXX = cumsum(PXX.*manyGW.*manyDelta);
    intPXY = cumsum(PXY.*manyGW.*manyDelta);
    intPYX = cumsum(PYX.*manyGW.*manyDelta);
    intPYY = cumsum(PYY.*manyGW.*manyDelta);
    
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
    subPXX(2:N,:) = intPXX(6:6:end-6,:);
    subPXY(2:N,:) = intPXY(6:6:end-6,:);
    subPYX(2:N,:) = intPYX(6:6:end-6,:);
    subPYY(2:N,:) = intPYY(6:6:end-6,:);
    
    GXX = intSXX(6:6:end,:)'-subSXX' + diag(2*(-2*deltaL/2.*log(deltaL/2)+3*deltaL/2));
    GXY = intSXY(6:6:end,:)'-subSXY';
    GYX = intSYX(6:6:end,:)'-subSYX';
    GYY = intSYY(6:6:end,:)'-subSYY' + diag(2*(-2*deltaL/2.*log(deltaL/2)+3*deltaL/2));
    TXXX = intQXXX(6:6:end,:)'-subQXXX';
    TXXY = intQXXY(6:6:end,:)'-subQXXY';
    TXYX = intQXYX(6:6:end,:)'-subQXYX';
    TXYY = intQXYY(6:6:end,:)'-subQXYY';
    TYXX = intQYXX(6:6:end,:)'-subQYXX';
    TYXY = intQYXY(6:6:end,:)'-subQYXY';
    TYYX = intQYYX(6:6:end,:)'-subQYYX';
    TYYY = intQYYY(6:6:end,:)'-subQYYY';
    DXX = intPXX(6:6:end,:)'-subPXX';
    DXY = intPXY(6:6:end,:)'-subPXY';
    DYX = intPYX(6:6:end,:)'-subPYX';
    DYY = intPYY(6:6:end,:)'-subPYY';
    
    R1 = repmat(r(1,:),N,1);
    R2 = repmat(r(2,:),N,1);
    
    A11 = TXXX.*R1+TXXY.*R2;
    A12 = TXYX.*R1+TXYY.*R2;
    A21 = TYXX.*R1+TYXY.*R2;
    A22 = TYYX.*R1+TYYY.*R2;
    D1 = DXX.*R1+DXY.*R2;
    D2 = DYX.*R1+DYY.*R2;
    
    %T = toc

end