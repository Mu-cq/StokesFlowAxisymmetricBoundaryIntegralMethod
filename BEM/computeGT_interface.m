%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,D1,D2,T1,T2,Iaxis] = computeGT_interface(x,y,walls,r)

    %tic

    %number of elements (than singularities)
    N = numel(x)-2;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    T1 = zeros(N);
    T2 = zeros(N);
    D1 = zeros(N);
    D2 = zeros(N);
    
    SXX = zeros(6*N,N);
    SXY = zeros(6*N,N);
    SYX = zeros(6*N,N);
    SYY = zeros(6*N,N);
    QXXX = zeros(6*N,N);
    QXXY = zeros(6*N,N);
    QXYX = zeros(6*N,N);
    QXYY = zeros(6*N,N);
    QYXX = zeros(6*N,N);
    QYXY = zeros(6*N,N);
    QYYX = zeros(6*N,N);
    QYYY = zeros(6*N,N);
    PXX = zeros(6*N,N);
    PXY = zeros(6*N,N);
    PYX = zeros(6*N,N);
    PYY = zeros(6*N,N);
    Iaxis = zeros(6*N,N);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = [(x(1:walls)+x(2:walls+1))/2 (x(walls+2:end-1)+x(walls+3:end))/2];
    Y0 = [(y(1:walls)+y(2:walls+1))/2 (y(walls+2:end-1)+y(walls+3:end))/2];
    
    %moltiplicate X0 and Y0
    XX(1:6:6*N-5) = X0;
    XX(2:6:6*N-4) = X0;
    XX(3:6:6*N-3) = X0;
    XX(4:6:6*N-2) = X0;
    XX(5:6:6*N-1) = X0;
    XX(6:6:6*N) = X0;
    YY(1:6:6*N-5) = Y0;
    YY(2:6:6*N-4) = Y0;
    YY(3:6:6*N-3) = Y0;
    YY(4:6:6*N-2) = Y0;
    YY(5:6:6*N-1) = Y0;
    YY(6:6:6*N) = Y0;
    
    % points where I perform gauss integration
    deltaX = [x(2:walls+1)-x(1:walls) x(walls+3:end)-x(walls+2:end-1)];
    deltaY = [y(2:walls+1)-y(1:walls) y(walls+3:end)-y(walls+2:end-1)];
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX(1:6:6*N-5) = GP(1)*deltaX/2;
    GPX(2:6:6*N-4) = GP(2)*deltaX/2;
    GPX(3:6:6*N-3) = GP(3)*deltaX/2;
    GPX(4:6:6*N-2) = GP(4)*deltaX/2;
    GPX(5:6:6*N-1) = GP(5)*deltaX/2;
    GPX(6:6:6*N) = GP(6)*deltaX/2;
    GPY(1:6:6*N-5) = GP(1)*deltaY/2;
    GPY(2:6:6*N-4) = GP(2)*deltaY/2;
    GPY(3:6:6*N-3) = GP(3)*deltaY/2;
    GPY(4:6:6*N-2) = GP(4)*deltaY/2;
    GPY(5:6:6*N-1) = GP(5)*deltaY/2;
    GPY(6:6:6*N) = GP(6)*deltaY/2;
    globalX = XX+GPX;
    globalY = YY+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    for k = 1:N
        for i = 1:6*N
            %not clear what is X0!!!
            [SXX(i,k),SXY(i,k),SYX(i,k),SYY(i,k),QXXX(i,k),QXXY(i,k),QXYX(i,k),QXYY(i,k),QYXX(i,k),QYXY(i,k),QYYX(i,k),QYYY(i,k),PXX(i,k),PXY(i,k),PYX(i,k),PYY(i,k) ,Iaxis(i,k)]...
                = sgf_ax_fs (2,globalX(i),globalY(i),X0(k),Y0(k));
            
        end
        %singularity treatment
        %display(k)
        SXX(1+6*(k-1):6*k,k) = SXX(1+6*(k-1):6*k,k)+2*log(sqrt((globalX(1+6*(k-1):6*k)-X0(k)).^2+(globalY(1+6*(k-1):6*k)-Y0(k)).^2)')-1;
        SYY(1+6*(k-1):6*k,k) = SYY(1+6*(k-1):6*k,k)+2*log(sqrt((globalX(1+6*(k-1):6*k)-X0(k)).^2+(globalY(1+6*(k-1):6*k)-Y0(k)).^2)')-1;
    end
    
    
    for i = 1:N
        for k = 1:N
            GXX(i,k) = SXX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            GXY(i,k) = SXY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            GYX(i,k) = SYX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            GYY(i,k) = SYY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            
            A11(i,k) = (QXXX(1+6*(k-1):6*k,i)'*r(1,k) + QXXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            A12(i,k) = (QXYX(1+6*(k-1):6*k,i)'*r(1,k) + QXYY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            A21(i,k) = (QYXX(1+6*(k-1):6*k,i)'*r(1,k) + QYXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            A22(i,k) = (QYYX(1+6*(k-1):6*k,i)'*r(1,k) + QYYY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            
            T1(i,k) = (QXXX(1+6*(k-1):6*k,i)'*r(1,k) + QXXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            T2(i,k) = (QYXX(1+6*(k-1):6*k,i)'*r(1,k) + QYXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            
            D1(i,k) = (PXX(1+6*(k-1):6*k,i)'*r(1,k) + PXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            D2(i,k) = (PYX(1+6*(k-1):6*k,i)'*r(1,k) + PYY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            
        end
        % singularity treatment
        GXX(i,i) = GXX(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
        GYY(i,i) = GYY(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
    end
    
    %T = toc

end