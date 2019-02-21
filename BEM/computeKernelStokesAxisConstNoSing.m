%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelStokesAxisConstNoSing(x,y,x0,y0,kernelFreeSpace,wallPos)

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x0';
    Y0 = y0';
    
    %compute normal vectors
    no = sqrt(diff(x).^2+diff(y).^2);
    nx = diff(y)./no;
    ny = -diff(x)./no;
    nnx = repmat(nx,6,1);
    nny = repmat(ny,6,1);
    nnx = repmat(nnx(:),1,N);
    nny = repmat(nny(:),1,N);
    
    %moltiplicate middle points coordinates
    xMiddle = (x(1:end-1)+x(2:end))/2;  yMiddle = (y(1:end-1)+y(2:end))/2;
    tempx = repmat(xMiddle,6,1);
    XX = reshape(tempx,1,6*elem);
    tempy = repmat(yMiddle,6,1);
    YY = reshape(tempy,1,6*elem);
    
    %lenght of the elements
    deltaX = diff(x);
    deltaY = diff(y);
    deltaL = sqrt(deltaX.^2+deltaY.^2);
    
    %every Gauss point
    GPX = repmat(GP,1,elem).*reshape((repmat(deltaX/2,6,1)),1,6*elem);
    GPY = repmat(GP,1,elem).*reshape((repmat(deltaY/2,6,1)),1,6*elem);
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
    X0matr = repmat(X0,6*elem,1);
    Y0matr = repmat(Y0,6*elem,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %compute gree function and the gradient
    if kernelFreeSpace==1     %free space
          [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
            sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
    elseif kernelFreeSpace==2     %wall bounded
        
          [SXX,SXY,SYX,SYY ...
            ,QXXX,QXXY,QXYX,QXYY ...
            ,QYXX,QYXY,QYYX,QYYY] = sgf_ax_w_vect (globalXmatr,globalYmatr,X0matr,Y0matr,wallPos);
          
    else
        error('Not implemented')
    end
    
    %mutuply kernnel of doublelayer times the normal
    Q11 = QXXX.*nnx + QXXY.*nny;
    Q12 = QXYX.*nnx + QXYY.*nny;
    Q21 = QYXX.*nnx + QYXY.*nny;
    Q22 = QYYX.*nnx + QYYY.*nny;
    
    temp = reshape(repmat(deltaL/2,6,1),1,6*elem);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',elem,N);
    
    %INTEGRATION
    intSXX = cumsum(SXX.*manyGW.*manyDelta);
    intSXY = cumsum(SXY.*manyGW.*manyDelta);
    intSYX = cumsum(SYX.*manyGW.*manyDelta);
    intSYY = cumsum(SYY.*manyGW.*manyDelta);
    intQ11 = cumsum(Q11.*manyGW.*manyDelta);
    intQ12 = cumsum(Q12.*manyGW.*manyDelta);
    intQ21 = cumsum(Q21.*manyGW.*manyDelta);
    intQ22 = cumsum(Q22.*manyGW.*manyDelta);
    
    subSXX(2:elem,:) = intSXX(6:6:end-6,:);
    subSXY(2:elem,:) = intSXY(6:6:end-6,:);
    subSYX(2:elem,:) = intSYX(6:6:end-6,:);
    subSYY(2:elem,:) = intSYY(6:6:end-6,:);
    subQ11(2:elem,:) = intQ11(6:6:end-6,:);
    subQ12(2:elem,:) = intQ12(6:6:end-6,:);
    subQ21(2:elem,:) = intQ21(6:6:end-6,:);
    subQ22(2:elem,:) = intQ22(6:6:end-6,:);
    
    GXX = intSXX(6:6:end,:)'-subSXX';
    GXY = intSXY(6:6:end,:)'-subSXY';
    GYX = intSYX(6:6:end,:)'-subSYX';
    GYY = intSYY(6:6:end,:)'-subSYY';
    A11 = intQ11(6:6:end,:)'-subQ11';
    A12 = intQ12(6:6:end,:)'-subQ12';
    A21 = intQ21(6:6:end,:)'-subQ21';
    A22 = intQ22(6:6:end,:)'-subQ22';

end