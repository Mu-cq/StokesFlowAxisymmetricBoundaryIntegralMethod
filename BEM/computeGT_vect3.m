%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_vect3(x,y)

    %tic

    %number of singularities
    N = numel(x)-1;
    
    %normal vector
    no = sqrt(diff(y(1:end)).^2+diff(x(1:end)).^2);
    r = [diff(y(1:end))./no; -diff(x(1:end))./no];
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = (x(1:end-1)+x(2:end))/2;
    Y0 = (y(1:end-1)+y(2:end))/2;
    
    %moltiplicate sing points coordinates
    tempx = repmat(X0,6,1);
    XX = reshape(tempx,1,6*N);
    tempy = repmat(Y0,6,1);
    YY = reshape(tempy,1,6*N);
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
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
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    %normal vecotr
    R1 = repmat(r(1,1:end),6,1);
    R2 = repmat(r(2,1:end),6,1);
    R1 = repmat(R1(:),1,N);
    R2 = repmat(R2(:),1,N);
            
    %kernel double layer
    T11 = QXXX.*R1 + QXXY.*R2;
    T12 = QXYX.*R1 + QXYY.*R2;
    T21 = QYXX.*R1 + QYXY.*R2;
    T22 = QYYX.*R1 + QYYY.*R2;
            
    for i = 1:N        
        %singularity treatment
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
    end

    temp = reshape(repmat(deltaL/2,6,1),1,6*N);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',N,N);
    
    %INTEGRATION
    intSXX = cumsum(SXX.*manyGW.*manyDelta);
    intSXY = cumsum(SXY.*manyGW.*manyDelta);
    intSYX = cumsum(SYX.*manyGW.*manyDelta);
    intSYY = cumsum(SYY.*manyGW.*manyDelta);
    intT11 = cumsum(T11.*manyGW.*manyDelta);
    intT12 = cumsum(T12.*manyGW.*manyDelta);
    intT21 = cumsum(T21.*manyGW.*manyDelta);
    intT22 = cumsum(T22.*manyGW.*manyDelta);
    
    subSXX(2:N,:) = intSXX(6:6:end-6,:);
    subSXY(2:N,:) = intSXY(6:6:end-6,:);
    subSYX(2:N,:) = intSYX(6:6:end-6,:);
    subSYY(2:N,:) = intSYY(6:6:end-6,:);
    subT11(2:N,:) = intT11(6:6:end-6,:);
    subT12(2:N,:) = intT12(6:6:end-6,:);
    subT21(2:N,:) = intT21(6:6:end-6,:);
    subT22(2:N,:) = intT22(6:6:end-6,:);
    
    GXX = intSXX(6:6:end,:)'-subSXX'+diag(2*(-2*deltaL/2.*log(deltaL/2)+3*deltaL/2));
    GXY = intSXY(6:6:end,:)'-subSXY';
    GYX = intSYX(6:6:end,:)'-subSYX';
    GYY = intSYY(6:6:end,:)'-subSYY'+diag(2*(-2*deltaL/2.*log(deltaL/2)+3*deltaL/2));
    A11 = intT11(6:6:end,:)'-subT11';
    A12 = intT12(6:6:end,:)'-subT12';
    A21 = intT21(6:6:end,:)'-subT21';
    A22 = intT22(6:6:end,:)'-subT22';
    
    %T = toc

end