%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY,PXX,PXY,PYX,PYY ,Iaxis] = computeGT(x,y)

    %tic

    %number of singularities
    N = numel(x)-1;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    TXXX = zeros(N);
    TXXY = zeros(N);
    TXYX = zeros(N);
    TXYY = zeros(N);
    TYXX = zeros(N);
    TYXY = zeros(N);
    TYYX = zeros(N);
    TYYY = zeros(N);
    
    SXX = zeros(numel(x)-1,6*(numel(x)-1));
    SXY = zeros(numel(x)-1,6*(numel(x)-1));
    SYX = zeros(numel(x)-1,6*(numel(x)-1));
    SYY = zeros(numel(x)-1,6*(numel(x)-1));
    QXXX = zeros(numel(x)-1,6*(numel(x)-1));
    QXXY = zeros(numel(x)-1,6*(numel(x)-1));
    QXYX = zeros(numel(x)-1,6*(numel(x)-1));
    QXYY = zeros(numel(x)-1,6*(numel(x)-1));
    QYXX = zeros(numel(x)-1,6*(numel(x)-1));
    QYXY = zeros(numel(x)-1,6*(numel(x)-1));
    QYYX = zeros(numel(x)-1,6*(numel(x)-1));
    QYYY = zeros(numel(x)-1,6*(numel(x)-1));
    PXX = zeros(numel(x)-1,6*(numel(x)-1));
    PXY = zeros(numel(x)-1,6*(numel(x)-1));
    PYX = zeros(numel(x)-1,6*(numel(x)-1));
    PYY = zeros(numel(x)-1,6*(numel(x)-1));
    Iaxis = zeros(numel(x)-1,6*(numel(x)-1));
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = (x(1:end-1)+x(2:end))/2;
    Y0 = (y(1:end-1)+y(2:end))/2;
    
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
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
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
    
    for i = 1:N
        for k = 1:6*N
            %not clear what is X0!!!
            [SXX(i,k),SXY(i,k),SYX(i,k),SYY(i,k),QXXX(i,k),QXXY(i,k),QXYX(i,k),QXYY(i,k),QYXX(i,k),QYXY(i,k),QYYX(i,k),QYYY(i,k),PXX(i,k),PXY(i,k),PYX(i,k),PYY(i,k) ,Iaxis(i,k)] = sgf_ax_fs (2,globalX(k),globalY(k),X0(i),Y0(i));
            
        end
        %singularity treatment
        SXX(i,1+6*(i-1):6*i) = SXX(i,1+6*(i-1):6*i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))-1;
        SYY(i,1+6*(i-1):6*i) = SYY(i,1+6*(i-1):6*i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))-1;
    end
    
    
    for i = 1:N
        for k = 1:N
            GXX(i,k) = SXX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GXY(i,k) = SXY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GYX(i,k) = SYX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GYY(i,k) = SYY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TXXX(i,k) = QXXX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TXXY(i,k) = QXXY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TXYX(i,k) = QXYX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TXYY(i,k) = QXYY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TYXX(i,k) = QYXX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TYXY(i,k) = QYXY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TYYX(i,k) = QYYX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            TYYY(i,k) = QYYY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
        end
        % singularity treatment
        GXX(i,i) = GXX(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
        GYY(i,i) = GYY(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
    end
    
    %T = toc

end