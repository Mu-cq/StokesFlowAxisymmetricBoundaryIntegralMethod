%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2,Iaxis] = computeGT_interface_lin(x,y,walls,r,inlet,outlet)

    %tic

    %total number of elements
    N = numel(x)-2;
    sing = N+1;
    
    GXX = zeros(N+1);
    GXY = zeros(N+1);
    GYX = zeros(N+1);
    GYY = zeros(N+1);
    A11 = zeros(N+1);
    A12 = zeros(N+1);
    A21 = zeros(N+1);
    A22 = zeros(N+1);
    T1 = zeros(N+1);
    T2 = zeros(N+1);
    D1 = zeros(N+1);
    D2 = zeros(N+1);
    
    SXX = zeros(6*N,N+1);
    SXY = zeros(6*N,N+1);
    SYX = zeros(6*N,N+1);
    SYY = zeros(6*N,N+1);
    QXXX = zeros(6*N,N+1);
    QXXY = zeros(6*N,N+1);
    QXYX = zeros(6*N,N+1);
    QXYY = zeros(6*N,N+1);
    QYXX = zeros(6*N,N+1);
    QYXY = zeros(6*N,N+1);
    QYYX = zeros(6*N,N+1);
    QYYY = zeros(6*N,N+1);
    PXX = zeros(6*N,N+1);
    PXY = zeros(6*N,N+1);
    PYX = zeros(6*N,N+1);
    PYY = zeros(6*N,N+1);
    Iaxis = zeros(6*N,N+1);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    phia = 1-(GP+1)/2;
    phib = (GP+1)/2;
    
    % point where the variable will be stored
    X0 = [(x(1:walls)+x(2:walls+1))/2 x(walls+2:end)];
    Y0 = [(y(1:walls)+y(2:walls+1))/2 y(walls+2:end)];
    
    %moltiplicate X0 and Y0
    XX(1:6:6*sing-5) = X0;
    XX(2:6:6*sing-4) = X0;
    XX(3:6:6*sing-3) = X0;
    XX(4:6:6*sing-2) = X0;
    XX(5:6:6*sing-1) = X0;
    XX(6:6:6*sing) = X0;
    YY(1:6:6*sing-5) = Y0;
    YY(2:6:6*sing-4) = Y0;
    YY(3:6:6*sing-3) = Y0;
    YY(4:6:6*sing-2) = Y0;
    YY(5:6:6*sing-1) = Y0;
    YY(6:6:6*sing) = Y0;
    
    % points where I perform gauss integration
    deltaX = [x(2:walls+1)-x(1:walls) x(walls+3:end)-x(walls+2:end-1)];
    deltaY = [y(2:walls+1)-y(1:walls) y(walls+3:end)-y(walls+2:end-1)];
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point (constant elements)
    GPX(1:6:6*N-5) = [GP(1)*deltaX(1:walls)/2 (GP(1)+1)*deltaX(walls+1:end)/2];
    GPX(2:6:6*N-4) = [GP(2)*deltaX(1:walls)/2 (GP(2)+1)*deltaX(walls+1:end)/2];
    GPX(3:6:6*N-3) = [GP(3)*deltaX(1:walls)/2 (GP(3)+1)*deltaX(walls+1:end)/2];
    GPX(4:6:6*N-2) = [GP(4)*deltaX(1:walls)/2 (GP(4)+1)*deltaX(walls+1:end)/2];
    GPX(5:6:6*N-1) = [GP(5)*deltaX(1:walls)/2 (GP(5)+1)*deltaX(walls+1:end)/2];
    GPX(6:6:6*N) = [GP(6)*deltaX(1:walls)/2 (GP(6)+1)*deltaX(walls+1:end)/2];
    GPY(1:6:6*N-5) = [GP(1)*deltaY(1:walls)/2 (GP(1)+1)*deltaY(walls+1:end)/2];
    GPY(2:6:6*N-4) = [GP(2)*deltaY(1:walls)/2 (GP(2)+1)*deltaY(walls+1:end)/2];
    GPY(3:6:6*N-3) = [GP(3)*deltaY(1:walls)/2 (GP(3)+1)*deltaY(walls+1:end)/2];
    GPY(4:6:6*N-2) = [GP(4)*deltaY(1:walls)/2 (GP(4)+1)*deltaY(walls+1:end)/2];
    GPY(5:6:6*N-1) = [GP(5)*deltaY(1:walls)/2 (GP(5)+1)*deltaY(walls+1:end)/2];
    GPY(6:6:6*N) = [GP(6)*deltaY(1:walls)/2 (GP(6)+1)*deltaY(walls+1:end)/2];
    
    globalX = XX(1:end-6)+GPX;
    globalY = YY(1:end-6)+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    for k = 1:N+1
        for i = 1:6*N
            %not clear what is X0!!!
            [SXX(i,k),SXY(i,k),SYX(i,k),SYY(i,k),QXXX(i,k),QXXY(i,k),QXYX(i,k),QXYY(i,k),QYXX(i,k),QYXY(i,k),QYYX(i,k),QYYY(i,k),PXX(i,k),PXY(i,k),PYX(i,k),PYY(i,k) ,Iaxis(i,k)]...
                = sgf_ax_fs (2,globalX(i),globalY(i),X0(k),Y0(k));
            
        end
        %singularity treatment
    end
    
    % singularity treatment const elem
    for k = 1:walls
        
        SXX(1+6*(k-1):6*k,k) = SXX(1+6*(k-1):6*k,k)+2*log(sqrt((globalX(1+6*(k-1):6*k)-X0(k)).^2+(globalY(1+6*(k-1):6*k)-Y0(k)).^2)')-1;
        SYY(1+6*(k-1):6*k,k) = SYY(1+6*(k-1):6*k,k)+2*log(sqrt((globalX(1+6*(k-1):6*k)-X0(k)).^2+(globalY(1+6*(k-1):6*k)-Y0(k)).^2)')-1;
        
    end
    
    % singularity treatment lin elem
     for i = walls+1:N
        
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
        SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
              
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))'-1;
        SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))'-1;
        
     end
    
     %add treatment for first and last node of the interface because they
     %are on the inlet and ouulet respectively
%      SXX(1+6*(inlet-1):6*inlet,walls+1) = SXX(1+6*(inlet-1):6*inlet,walls+1) + 2*log(sqrt((globalX(1+6*(inlet-1):6*inlet)-X0(walls+1)).^2+(globalY(1+6*(inlet-1):6*inlet)-Y0(walls+1)).^2)')-1;
%      SYY(1+6*(inlet-1):6*inlet,walls+1) = SYY(1+6*(inlet-1):6*inlet,walls+1) + 2*log(sqrt((globalX(1+6*(inlet-1):6*inlet)-X0(walls+1)).^2+(globalY(1+6*(inlet-1):6*inlet)-Y0(walls+1)).^2)')-1;
%      SXX(1+6*(inlet):6*(inlet+1),walls+1) = SXX(1+6*(inlet):6*(inlet+1),walls+1) + 2*log(sqrt((globalX(1+6*(inlet):6*(inlet+1))-X0(walls+1)).^2+(globalY(1+6*(inlet):6*(inlet+1))-Y0(walls+1)).^2)')-1;
%      SYY(1+6*(inlet):6*(inlet+1),walls+1) = SYY(1+6*(inlet):6*(inlet+1),walls+1) + 2*log(sqrt((globalX(1+6*(inlet):6*(inlet+1))-X0(walls+1)).^2+(globalY(1+6*(inlet):6*(inlet+1))-Y0(walls+1)).^2)')-1;
%      
%      SXX(1+6*(outlet-1):6*outlet,end) = SXX(1+6*(outlet-1):6*outlet,end) + 2*log(sqrt((globalX(1+6*(outlet-1):6*outlet)-X0(end)).^2+(globalY(1+6*(outlet-1):6*outlet)-Y0(end)).^2)')-1;
%      SYY(1+6*(outlet-1):6*outlet,end) = SYY(1+6*(outlet-1):6*outlet,end) + 2*log(sqrt((globalX(1+6*(outlet-1):6*outlet)-X0(end)).^2+(globalY(1+6*(outlet-1):6*outlet)-Y0(end)).^2)')-1;
%      SXX(1+6*(outlet):6*(outlet+1),end) = SXX(1+6*(outlet):6*(outlet+1),end) + 2*log(sqrt((globalX(1+6*(outlet):6*(outlet+1))-X0(end)).^2+(globalY(1+6*(outlet):6*(outlet+1))-Y0(end)).^2)')-1;
%      SYY(1+6*(outlet):6*(outlet+1),end) = SYY(1+6*(outlet):6*(outlet+1),end) + 2*log(sqrt((globalX(1+6*(outlet):6*(outlet+1))-X0(end)).^2+(globalY(1+6*(outlet):6*(outlet+1))-Y0(end)).^2)')-1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT ELEMENTS INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:N+1
        for k = 1:walls
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
%         GXX(i,i) = GXX(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
%         GYY(i,i) = GYY(i,i)+2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
    end
    
    % singularity treatment
    for i = 1:walls
        GXX(i,i) = GXX(i,i) + 2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
        GYY(i,i) = GYY(i,i) + 2*(-2*deltaL(i)/2*log(deltaL(i)/2)+3*deltaL(i)/2);
    end
    
    %add treatment for first and last node of the interface because they
     %are on the inlet and ouulet respectively
%     GXX(walls+1,inlet) = GXX(walls+1,inlet) + 2*(-2*deltaL(inlet)/2*log(deltaL(inlet)/2)+3*deltaL(inlet)/2);
%     GYY(walls+1,inlet) = GYY(walls+1,inlet) + 2*(-2*deltaL(inlet)/2*log(deltaL(inlet)/2)+3*deltaL(inlet)/2);
%     GXX(walls+1,inlet+1) = GXX(walls+1,inlet+1) + 2*(-2*deltaL(inlet+1)/2*log(deltaL(inlet+1)/2)+3*deltaL(inlet+1)/2);
%     GYY(walls+1,inlet+1) = GYY(walls+1,inlet+1) + 2*(-2*deltaL(inlet+1)/2*log(deltaL(inlet+1)/2)+3*deltaL(inlet+1)/2);
%     
%     GXX(end,outlet) = GXX(end,outlet) + 2*(-2*deltaL(outlet)/2*log(deltaL(outlet)/2)+3*deltaL(outlet)/2);
%     GYY(end,outlet) = GYY(end,outlet) + 2*(-2*deltaL(outlet)/2*log(deltaL(outlet)/2)+3*deltaL(outlet)/2);
%     GXX(end,outlet+1) = GXX(end,outlet+1) + 2*(-2*deltaL(outlet+1)/2*log(deltaL(outlet+1)/2)+3*deltaL(outlet+1)/2);
%     GYY(end,outlet+1) = GYY(end,outlet+1) + 2*(-2*deltaL(outlet+1)/2*log(deltaL(outlet+1)/2)+3*deltaL(outlet+1)/2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR ELEMENTS INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     for i = 1:N
        for k = walls+1:N
            GXX(i,k) = GXX(i,k) + SXX(1+6*(k-1):6*k,i)'.*phia*GW'*deltaL(k)/2;
            GXX(i,k+1) = SXX(1+6*(k-1):6*k,i)'.*phib*GW'*deltaL(k)/2;
            
            GXY(i,k) = GXY(i,k) + SXY(1+6*(k-1):6*k,i)'.*phia*GW'*deltaL(k)/2;
            GXY(i,k+1) = SXY(1+6*(k-1):6*k,i)'.*phib*GW'*deltaL(k)/2;
            
            GYX(i,k) = GYX(i,k) + SYX(1+6*(k-1):6*k,i)'.*phia*GW'*deltaL(k)/2;
            GYX(i,k+1) = SYX(1+6*(k-1):6*k,i)'.*phib*GW'*deltaL(k)/2;
            
            GYY(i,k) = GYY(i,k) + SYY(1+6*(k-1):6*k,i)'.*phia*GW'*deltaL(k)/2;
            GYY(i,k+1) = SYY(1+6*(k-1):6*k,i)'.*phib*GW'*deltaL(k)/2;
            
            A11(i,k) = A11(i,k) + (QXXX(1+6*(k-1):6*k,i)'*r(1,k) + QXXY(1+6*(k-1):6*k,i)'*r(2,k)).*phia*GW'*deltaL(k)/2;
            A11(i,k+1) = (QXXX(1+6*(k-1):6*k,i)'*r(1,k) + QXXY(1+6*(k-1):6*k,i)'*r(2,k)).*phib*GW'*deltaL(k)/2;
            
            A12(i,k) = A12(i,k) + (QXYX(1+6*(k-1):6*k,i)'*r(1,k) + QXYY(1+6*(k-1):6*k,i)'*r(2,k)).*phia*GW'*deltaL(k)/2;
            A12(i,k+1) = (QXYX(1+6*(k-1):6*k,i)'*r(1,k) + QXYY(1+6*(k-1):6*k,i)'*r(2,k)).*phib*GW'*deltaL(k)/2;
            
            A21(i,k) = A21(i,k) + (QYXX(1+6*(k-1):6*k,i)'*r(1,k) + QYXY(1+6*(k-1):6*k,i)'*r(2,k)).*phia*GW'*deltaL(k)/2;
            A21(i,k+1) = (QYXX(1+6*(k-1):6*k,i)'*r(1,k) + QYXY(1+6*(k-1):6*k,i)'*r(2,k)).*phib*GW'*deltaL(k)/2;
            
            A22(i,k) = A22(i,k) + (QYYX(1+6*(k-1):6*k,i)'*r(1,k) + QYYY(1+6*(k-1):6*k,i)'*r(2,k)).*phia*GW'*deltaL(k)/2;
            A22(i,k+1) = (QYYX(1+6*(k-1):6*k,i)'*r(1,k) + QYYY(1+6*(k-1):6*k,i)'*r(2,k)).*phib*GW'*deltaL(k)/2;
            
            T1(i,k) = (QXXX(1+6*(k-1):6*k,i)'*r(1,k) + QXXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            T2(i,k) = (QYXX(1+6*(k-1):6*k,i)'*r(1,k) + QYXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            
            D1(i,k) = (PXX(1+6*(k-1):6*k,i)'*r(1,k) + PXY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
            D2(i,k) = (PYX(1+6*(k-1):6*k,i)'*r(1,k) + PYY(1+6*(k-1):6*k,i)'*r(2,k))*GW'*deltaL(k)/2;
        end
    end
    
    %singularity treatment for linear element REPLACE WITH VECTORIAL
    %PRODUCT
    for i = walls+1:N
        
       GXX(i,i) = GXX(i,i) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
       GXX(i,i+1) = GXX(i,i+1) - deltaL(i)*log(deltaL(i))+deltaL(i);
       GXX(i+1,i) = GXX(i+1,i) - deltaL(i)*log(deltaL(i))+deltaL(i);
       GXX(i+1,i+1) = GXX(i+1,i+1) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
       
       GYY(i,i) = GYY(i,i) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
       GYY(i,i+1) = GYY(i,i+1) - deltaL(i)*log(deltaL(i))+deltaL(i);
       GYY(i+1,i) = GYY(i+1,i) - deltaL(i)*log(deltaL(i))+deltaL(i);
       GYY(i+1,i+1) = GYY(i+1,i+1) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
       
    end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %T = toc

end