%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2,F11,F12,F21,F22,Iaxis] = computeGT_linear_more3(x,y,phia,phib,r)

    %tic

    %number of nodes
    N = numel(x);
    %number of element
    el = N-1;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    T1 = zeros(N,N-1);
    T2 = zeros(N,N-1);
    D1 = zeros(N,N-1);
    D2 = zeros(N,N-1);
    F11 = zeros(N,N-1);
    F12 = zeros(N,N-1);
    F21 = zeros(N,N-1);
    F22 = zeros(N,N-1);
    
    SXX = zeros(N-1,6*(N-1));
    SXY = zeros(N-1,6*(N-1));
    SYX = zeros(N-1,6*(N-1));
    SYY = zeros(N-1,6*(N-1));
    QXXX = zeros(N-1,6*(N-1));
    QXXY = zeros(N-1,6*(N-1));
    QXYX = zeros(N-1,6*(N-1));
    QXYY = zeros(N-1,6*(N-1));
    QYXX = zeros(N-1,6*(N-1));
    QYXY = zeros(N-1,6*(N-1));
    QYYX = zeros(N-1,6*(N-1));
    QYYY = zeros(N-1,6*(N-1));
    PXX = zeros(N-1,6*(N-1));
    PXY = zeros(N-1,6*(N-1));
    PYX = zeros(N-1,6*(N-1));
    PYY = zeros(N-1,6*(N-1));
    Iaxis = zeros(N-1,6*(N-1));
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x;
    Y0 = y;
    
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
    GP = GP+1;
    GPX(1:6:6*el-5) = GP(1)*deltaX/2;
    GPX(2:6:6*el-4) = GP(2)*deltaX/2;
    GPX(3:6:6*el-3) = GP(3)*deltaX/2;
    GPX(4:6:6*el-2) = GP(4)*deltaX/2;
    GPX(5:6:6*el-1) = GP(5)*deltaX/2;
    GPX(6:6:6*el) = GP(6)*deltaX/2;
    GPY(1:6:6*el-5) = GP(1)*deltaY/2;
    GPY(2:6:6*el-4) = GP(2)*deltaY/2;
    GPY(3:6:6*el-3) = GP(3)*deltaY/2;
    GPY(4:6:6*el-2) = GP(4)*deltaY/2;
    GPY(5:6:6*el-1) = GP(5)*deltaY/2;
    GPY(6:6:6*el) = GP(6)*deltaY/2;
    globalX = XX(1:end-6)+GPX;
    globalY = YY(1:end-6)+GPY;
    
%     figure
%     plot(globalX,phia,globalX,phib)
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    for i = 1:N
        for k = 1:6*(N-1)
            %not clear what is X0!!!
            [SXX(i,k),SXY(i,k),SYX(i,k),SYY(i,k),QXXX(i,k),QXXY(i,k),QXYX(i,k),QXYY(i,k),QYXX(i,k),QYXY(i,k),QYYX(i,k),QYYY(i,k),PXX(i,k),PXY(i,k),PYX(i,k),PYY(i,k) ,Iaxis(i,k)] = sgf_ax_fs (2,globalX(k),globalY(k),X0(i),Y0(i));
        end
    end
    
    %singularity treatment for linear element REPLACE WITH VECTORIAL
    %PRODUCT %%%%%
%     for i = 2:N-2
%        SXX(i,1+6*(i-1):6*i) = SXX(i,1+6*(i-1):6*i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))-1;
%        SXX(i+1,1+6*(i-1):6*i) = SXX(i+1,1+6*(i-1):6*i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))-1;
%               
%        SYY(i,1+6*(i-1):6*i) = SYY(i,1+6*(i-1):6*i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))-1;
%        SYY(i+1,1+6*(i-1):6*i) = SYY(i+1,1+6*(i-1):6*i) + 2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i+1)).^2+(globalY(1+6*(i-1):6*i)-Y0(i+1)).^2))-1;
%     end
    
%     SXX(2,1:6) = SXX(2,1:6) + 2*log(sqrt((globalX(1:6)-X0(2)).^2+(globalY(1:6)-Y0(2)).^2))-1;
%     SXX(N-1,1+6*(N-2):6*(N-1)) = SXX(N-1,1+6*(N-2):6*(N-1)) + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))-1;
%     
%     SYY(2,1:6) = SYY(2,1:6) + 2*log(sqrt((globalX(1:6)-X0(2)).^2+(globalY(1:6)-Y0(2)).^2))-1;
%     SYY(N-1,1+6*(N-2):6*(N-1)) = SYY(N-1,1+6*(N-2):6*(N-1)) + 2*log(sqrt((globalX(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(globalY(1+6*(N-2):6*(N-1))-Y0(N-1)).^2))-1;
    
    %vectorial products for sing treatment SEEMS SLOWER!!!
%     SXX(1:end-1,1:6:end-5) = SXX(1:end-1,1:6:end-5) + diag(2*log(sqrt((globalX(1:6:end-5)-XX(1:6:end-11)).^2+(globalY(1:6:end-5)-YY(1:6:end-11)).^2))-1);
%     SXX(1:end-1,2:6:end-4) = SXX(1:end-1,2:6:end-4) + diag(2*log(sqrt((globalX(2:6:end-4)-XX(2:6:end-10)).^2+(globalY(2:6:end-4)-YY(2:6:end-10)).^2))-1);
%     SXX(1:end-1,3:6:end-3) = SXX(1:end-1,3:6:end-3) + diag(2*log(sqrt((globalX(3:6:end-3)-XX(3:6:end-9)).^2+(globalY(3:6:end-3)-YY(3:6:end-9)).^2))-1);
%     SXX(1:end-1,4:6:end-2) = SXX(1:end-1,4:6:end-2) + diag(2*log(sqrt((globalX(4:6:end-2)-XX(4:6:end-8)).^2+(globalY(4:6:end-2)-YY(4:6:end-8)).^2))-1);
%     SXX(1:end-1,5:6:end-1) = SXX(1:end-1,5:6:end-1) + diag(2*log(sqrt((globalX(5:6:end-1)-XX(5:6:end-7)).^2+(globalY(5:6:end-1)-YY(5:6:end-7)).^2))-1);
%     SXX(1:end-1,6:6:end) = SXX(1:end-1,6:6:end) + diag(2*log(sqrt((globalX(6:6:end)-XX(6:6:end-6)).^2+(globalY(6:6:end)-YY(6:6:end-6)).^2))-1);
%     
%     SXX(2:end,1:6:end-5) = SXX(2:end,1:6:end-5) + diag(2*log(sqrt((globalX(1:6:end-5)-XX(7:6:end-5)).^2+(globalY(1:6:end-5)-YY(7:6:end-5)).^2))-1);
%     SXX(2:end,2:6:end-4) = SXX(2:end,2:6:end-4) + diag(2*log(sqrt((globalX(2:6:end-4)-XX(8:6:end-4)).^2+(globalY(2:6:end-4)-YY(8:6:end-4)).^2))-1);
%     SXX(2:end,3:6:end-3) = SXX(2:end,3:6:end-3) + diag(2*log(sqrt((globalX(3:6:end-3)-XX(9:6:end-3)).^2+(globalY(3:6:end-3)-YY(9:6:end-3)).^2))-1);
%     SXX(2:end,4:6:end-2) = SXX(2:end,4:6:end-2) + diag(2*log(sqrt((globalX(4:6:end-2)-XX(10:6:end-2)).^2+(globalY(4:6:end-2)-YY(10:6:end-2)).^2))-1);
%     SXX(2:end,5:6:end-1) = SXX(2:end,5:6:end-1) + diag(2*log(sqrt((globalX(5:6:end-1)-XX(11:6:end-1)).^2+(globalY(5:6:end-1)-YY(11:6:end-1)).^2))-1);
%     SXX(2:end,6:6:end) = SXX(2:end,6:6:end) + diag(2*log(sqrt((globalX(6:6:end)-XX(12:6:end)).^2+(globalY(6:6:end)-YY(12:6:end)).^2))-1);
%     
%     SYY(1:end-1,1:6:end-5) = SYY(1:end-1,1:6:end-5) + diag(2*log(sqrt((globalX(1:6:end-5)-XX(1:6:end-11)).^2+(globalY(1:6:end-5)-YY(1:6:end-11)).^2))-1);
%     SYY(1:end-1,2:6:end-4) = SYY(1:end-1,2:6:end-4) + diag(2*log(sqrt((globalX(2:6:end-4)-XX(2:6:end-10)).^2+(globalY(2:6:end-4)-YY(2:6:end-10)).^2))-1);
%     SYY(1:end-1,3:6:end-3) = SYY(1:end-1,3:6:end-3) + diag(2*log(sqrt((globalX(3:6:end-3)-XX(3:6:end-9)).^2+(globalY(3:6:end-3)-YY(3:6:end-9)).^2))-1);
%     SYY(1:end-1,4:6:end-2) = SYY(1:end-1,4:6:end-2) + diag(2*log(sqrt((globalX(4:6:end-2)-XX(4:6:end-8)).^2+(globalY(4:6:end-2)-YY(4:6:end-8)).^2))-1);
%     SYY(1:end-1,5:6:end-1) = SYY(1:end-1,5:6:end-1) + diag(2*log(sqrt((globalX(5:6:end-1)-XX(5:6:end-7)).^2+(globalY(5:6:end-1)-YY(5:6:end-7)).^2))-1);
%     SYY(1:end-1,6:6:end) = SYY(1:end-1,6:6:end) + diag(2*log(sqrt((globalX(6:6:end)-XX(6:6:end-6)).^2+(globalY(6:6:end)-YY(6:6:end-6)).^2))-1);
%     
%     SYY(2:end,1:6:end-5) = SYY(2:end,1:6:end-5) + diag(2*log(sqrt((globalX(1:6:end-5)-XX(7:6:end-5)).^2+(globalY(1:6:end-5)-YY(7:6:end-5)).^2))-1);
%     SYY(2:end,2:6:end-4) = SYY(2:end,2:6:end-4) + diag(2*log(sqrt((globalX(2:6:end-4)-XX(8:6:end-4)).^2+(globalY(2:6:end-4)-YY(8:6:end-4)).^2))-1);
%     SYY(2:end,3:6:end-3) = SYY(2:end,3:6:end-3) + diag(2*log(sqrt((globalX(3:6:end-3)-XX(9:6:end-3)).^2+(globalY(3:6:end-3)-YY(9:6:end-3)).^2))-1);
%     SYY(2:end,4:6:end-2) = SYY(2:end,4:6:end-2) + diag(2*log(sqrt((globalX(4:6:end-2)-XX(10:6:end-2)).^2+(globalY(4:6:end-2)-YY(10:6:end-2)).^2))-1);
%     SYY(2:end,5:6:end-1) = SYY(2:end,5:6:end-1) + diag(2*log(sqrt((globalX(5:6:end-1)-XX(11:6:end-1)).^2+(globalY(5:6:end-1)-YY(11:6:end-1)).^2))-1);
%     SYY(2:end,6:6:end) = SYY(2:end,6:6:end) + diag(2*log(sqrt((globalX(6:6:end)-XX(12:6:end)).^2+(globalY(6:6:end)-YY(12:6:end)).^2))-1);
    
    %integration for stresses and velocities, linear elements
    for i = 1:N
        for k = 1:N-1
            GXX(i,k) = GXX(i,k) + SXX(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GXX(i,k+1) = SXX(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            GXY(i,k) = GXY(i,k) + SXY(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GXY(i,k+1) = SXY(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            GYX(i,k) = GYX(i,k) + SYX(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GYX(i,k+1) = SYX(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            GYY(i,k) = GYY(i,k) + SYY(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            GYY(i,k+1) = SYY(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            A11(i,k) = A11(i,k) + (QXXX(i,1+6*(k-1):6*k)*r(1,k) + QXXY(i,1+6*(k-1):6*k)*r(2,k)).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            A11(i,k+1) = (QXXX(i,1+6*(k-1):6*k)*r(1,k) + QXXY(i,1+6*(k-1):6*k)*r(2,k)).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            A12(i,k) = A12(i,k) + (QXYX(i,1+6*(k-1):6*k)*r(1,k) + QXYY(i,1+6*(k-1):6*k)*r(2,k)).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            A12(i,k+1) = (QXYX(i,1+6*(k-1):6*k)*r(1,k) + QXYY(i,1+6*(k-1):6*k)*r(2,k)).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            A21(i,k) = A21(i,k) + (QYXX(i,1+6*(k-1):6*k)*r(1,k) + QYXY(i,1+6*(k-1):6*k)*r(2,k)).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            A21(i,k+1) = (QYXX(i,1+6*(k-1):6*k)*r(1,k) + QYXY(i,1+6*(k-1):6*k)*r(2,k)).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            A22(i,k) = A22(i,k) + (QYYX(i,1+6*(k-1):6*k)*r(1,k) + QYYY(i,1+6*(k-1):6*k)*r(2,k)).*phia(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            A22(i,k+1) = (QYYX(i,1+6*(k-1):6*k)*r(1,k) + QYYY(i,1+6*(k-1):6*k)*r(2,k)).*phib(1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            
            T1(i,k) = (QXXX(i,1+6*(k-1):6*k)*r(1,k) + QXXY(i,1+6*(k-1):6*k)*r(2,k))*GW'*deltaL(k)/2;
            T2(i,k) = (QYXX(i,1+6*(k-1):6*k)*r(1,k) + QYXY(i,1+6*(k-1):6*k)*r(2,k))*GW'*deltaL(k)/2;
            
            D1(i,k) = (PXX(i,1+6*(k-1):6*k)*r(1,k) + PXY(i,1+6*(k-1):6*k)*r(2,k))*GW'*deltaL(k)/2;
            D2(i,k) = (PYX(i,1+6*(k-1):6*k)*r(1,k) + PYY(i,1+6*(k-1):6*k)*r(2,k))*GW'*deltaL(k)/2;
            
            F11(i,k) = SXX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            F12(i,k) = SXY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            F21(i,k) = SYX(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
            F22(i,k) = SYY(i,1+6*(k-1):6*k)*GW'*deltaL(k)/2;
        end
    end
    
    %singularity treatment for linear element REPLACE WITH VECTORIAL
    %PRODUCT
%     for i = 2:N-2
%         
%        GXX(i,i) = GXX(i,i) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
%        GXX(i,i+1) = GXX(i,i+1) - deltaL(i)*log(deltaL(i))+deltaL(i);
%        GXX(i+1,i) = GXX(i+1,i) - deltaL(i)*log(deltaL(i))+deltaL(i);
%        GXX(i+1,i+1) = GXX(i+1,i+1) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
%        
%        GYY(i,i) = GYY(i,i) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
%        GYY(i,i+1) = GYY(i,i+1) - deltaL(i)*log(deltaL(i))+deltaL(i);
%        GYY(i+1,i) = GYY(i+1,i) - deltaL(i)*log(deltaL(i))+deltaL(i);
%        GYY(i+1,i+1) = GYY(i+1,i+1) - deltaL(i)*log(deltaL(i))+2*deltaL(i);
%        
%     end
%     
%     GXX(2,1) = GXX(2,1) - deltaL(1)*log(deltaL(1))+deltaL(1);
%     GXX(2,2) = GXX(2,2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
%     GXX(N-1,N-1) = GXX(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
%     GXX(N-1,N) = GXX(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
%     
%     GYY(2,1) = GYY(2,1) - deltaL(1)*log(deltaL(1))+deltaL(1);
%     GYY(2,2) = GYY(2,2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
%     GYY(N-1,N-1) = GYY(N-1,N-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
%     GYY(N-1,N) = GYY(N-1,N) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);

    %vectorial products for sing treatment SEEMS SLOWER!!!
%     GXX(1:end-1,1:end-1) = GXX(1:end-1,1:end-1) + diag(-deltaL.*log(deltaL)+2*deltaL);
%     GXX(1:end-1,2:end) = GXX(1:end-1,2:end) + diag(-deltaL.*log(deltaL)+deltaL);
%     GXX(2:end,1:end-1) = GXX(2:end,1:end-1) + diag(-deltaL.*log(deltaL)+deltaL);
%     GXX(2:end,2:end) = GXX(2:end,2:end) + diag(-deltaL.*log(deltaL)+2*deltaL);
%     
%     GYY(1:end-1,1:end-1) = GYY(1:end-1,1:end-1) + diag(-deltaL.*log(deltaL)+2*deltaL);
%     GYY(1:end-1,2:end) = GYY(1:end-1,2:end) + diag(-deltaL.*log(deltaL)+deltaL);
%     GYY(2:end,1:end-1) = GYY(2:end,1:end-1) + diag(-deltaL.*log(deltaL)+deltaL);
%     GYY(2:end,2:end) = GYY(2:end,2:end) + diag(-deltaL.*log(deltaL)+2*deltaL);
    
    %T = toc

end