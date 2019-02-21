%compute Green's functions and associated stress tensor everywhere is
%needed for Brinkmann

function [GXX,GXY,GYY,TXXX,TXXY,TXYY,TYYY] = computeGT_2DStokes_vlado(x,y,wei)

    %tic

    %number of singularities
    N = numel(x)-1;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    %GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %weigths matrix
    %wei(1:,1:N) = ;
    
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
    [SXX,SXY,SYY,QXXX,QXXY,QXYY,QYYY] =...
               GT_2DStokes (globalXmatr,globalYmatr,X0matr,Y0matr);
            
%     xgauss = (repmat(X0',1,6)+repmat(deltaX',1,6))';
%     ygauss = (repmat(Y0',1,6)+repmat(deltaY',1,6))';
%     xgauss = xgauss(:)';
%     ygauss = ygauss(:)';
%     
%     [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY]...
%         = GT_2DStokes_cfunction(xgauss,ygauss,X0,Y0);
%     
%     SXX = reshape(SXX,numel(xgauss),numel(X0));
%     SXY = reshape(SXY,numel(xgauss),numel(X0));
%     SYX = reshape(SYX,numel(xgauss),numel(X0));
%     SYY = reshape(SYY,numel(xgauss),numel(X0));
%     QXXX = reshape(QXXX,numel(xgauss),numel(X0));
%     QXXY = reshape(QXXY,numel(xgauss),numel(X0));
%     QXYX = reshape(QXYX,numel(xgauss),numel(X0));
%     QXYY = reshape(QXYY,numel(xgauss),numel(X0));
%     QYXX = reshape(QYXX,numel(xgauss),numel(X0));
%     QYXY = reshape(QYXY,numel(xgauss),numel(X0));
%     QYYX = reshape(QYYX,numel(xgauss),numel(X0));
%     QYYY = reshape(QYYY,numel(xgauss),numel(X0));
    
            
    for i = 1:N        
        %singularity treatment
        %SINGLE LAYER
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i)+log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))';
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i)+log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))';
        
        %DOUBLE LAYER
        QXXX(1+6*(i-1):6*i,i) = zeros(6,1);
        QXXY(1+6*(i-1):6*i,i) = zeros(6,1);
        QXYY(1+6*(i-1):6*i,i) = zeros(6,1);
        QYYY(1+6*(i-1):6*i,i) = zeros(6,1);
    end

    temp = reshape(repmat(deltaL/2,6,1),1,6*N);
    manyDelta = repmat(temp',1,N);
    %manyGW = repmat(GW',N,N);
    
    %INTEGRATION    
    GXX = wei*(SXX.*manyDelta) + diag(-deltaL.*log(deltaL/2)+deltaL);
    GXY = wei*(SXY.*manyDelta);
    GYY = wei*(SYY.*manyDelta) + diag(-deltaL.*log(deltaL/2)+deltaL);
    TXXX = wei*(QXXX.*manyDelta);
    TXXY = wei*(QXXY.*manyDelta);
    TXYY = wei*(QXYY.*manyDelta);
    TYYY = wei*(QYYY.*manyDelta);
    
    %T = toc

end