%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY,PXX,PXY,PYX,PYY ,Iaxis] = computeGT_vect(x,y)

    tic

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
        
            %not clear what is X0!!!
            [SXX(:,i),SXY(:,i),SYX(:,i),SYY(:,i),QXXX(:,i),QXXY(:,i),QXYX(:,i),QXYY(:,i),QYXX(:,i),QYXY(:,i),QYYX(:,i),QYYY(:,i),PXX(:,i),PXY(:,i),PYX(:,i),PYY(:,i) ,Iaxis(:,i)] =...
                sgf_ax_fs_vect (2,globalX',globalY',X0(i)*ones(6*N,1),Y0(i)*ones(6*N,1));
            
        
        %singularity treatment
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
    end
    
    
    for i = 1:N
        for k = 1:N
            GXX(i,k) = SXX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            GXY(i,k) = SXY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            GYX(i,k) = SYX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            GYY(i,k) = SYY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TXXX(i,k) = QXXX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TXXY(i,k) = QXXY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TXYX(i,k) = QXYX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TXYY(i,k) = QXYY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TYXX(i,k) = QYXX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TYXY(i,k) = QYXY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TYYX(i,k) = QYYX(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
            TYYY(i,k) = QYYY(1+6*(k-1):6*k,i)'*GW'*deltaL(k)/2;
        end
        % singularity treatment
        GXX(i,i) = GXX(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
        GYY(i,i) = GYY(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
    end
    
    T = toc

end