%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2,Iaxis] = computeGT_spline2(x,y,phia,phib)

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
    
    SXX = zeros(el,6*el);
    SXY = zeros(el,6*el);
    SYX = zeros(el,6*el);
    SYY = zeros(el,6*el);
    QXXX = zeros(el,6*el);
    QXXY = zeros(el,6*el);
    QXYX = zeros(el,6*el);
    QXYY = zeros(el,6*el);
    QYXX = zeros(el,6*el);
    QYXY = zeros(el,6*el);
    QYYX = zeros(el,6*el);
    QYYY = zeros(el,6*el);
    PXX = zeros(el,6*el);
    PXY = zeros(el,6*el);
    PYX = zeros(el,6*el);
    PYY = zeros(el,6*el);
    Iaxis = zeros(el,6*el);
    
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (x, y);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    aax(1:6:6*el-5) = ax;
    aax(2:6:6*el-4) = ax;
    aax(3:6:6*el-3) = ax;
    aax(4:6:6*el-2) = ax;
    aax(5:6:6*el-1) = ax;
    aax(6:6:6*el) = ax;
    bbx(1:6:6*el-5) = bx;
    bbx(2:6:6*el-4) = bx;
    bbx(3:6:6*el-3) = bx;
    bbx(4:6:6*el-2) = bx;
    bbx(5:6:6*el-1) = bx;
    bbx(6:6:6*el) = bx;
    ccx(1:6:6*el-5) = cx;
    ccx(2:6:6*el-4) = cx;
    ccx(3:6:6*el-3) = cx;
    ccx(4:6:6*el-2) = cx;
    ccx(5:6:6*el-1) = cx;
    ccx(6:6:6*el) = cx;
    ddx(1:6:6*el-5) = dx;
    ddx(2:6:6*el-4) = dx;
    ddx(3:6:6*el-3) = dx;
    ddx(4:6:6*el-2) = dx;
    ddx(5:6:6*el-1) = dx;
    ddx(6:6:6*el) = dx;
    aay(1:6:6*el-5) = ay;
    aay(2:6:6*el-4) = ay;
    aay(3:6:6*el-3) = ay;
    aay(4:6:6*el-2) = ay;
    aay(5:6:6*el-1) = ay;
    aay(6:6:6*el) = ay;
    bby(1:6:6*el-5) = by;
    bby(2:6:6*el-4) = by;
    bby(3:6:6*el-3) = by;
    bby(4:6:6*el-2) = by;
    bby(5:6:6*el-1) = by;
    bby(6:6:6*el) = by;
    ccy(1:6:6*el-5) = cy;
    ccy(2:6:6*el-4) = cy;
    ccy(3:6:6*el-3) = cy;
    ccy(4:6:6*el-2) = cy;
    ccy(5:6:6*el-1) = cy;
    ccy(6:6:6*el) = cy;
    ddy(1:6:6*el-5) = dy;
    ddy(2:6:6*el-4) = dy;
    ddy(3:6:6*el-3) = dy;
    ddy(4:6:6*el-2) = dy;
    ddy(5:6:6*el-1) = dy;
    ddy(6:6:6*el) = dy;
    
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;
    GW = GW/2;

    % point where the variable will be stored
    X0 = x;
    Y0 = y;
    
    % points where I perform gauss integration
    GPglobal(1:6:6*el-5) = GPt(1);
    GPglobal(2:6:6*el-4) = GPt(2);
    GPglobal(3:6:6*el-3) = GPt(3);
    GPglobal(4:6:6*el-2) = GPt(4);
    GPglobal(5:6:6*el-1) = GPt(5);
    GPglobal(6:6:6*el) = GPt(6);
    
    %global coordianted retrived on the splines
    eta = aax+bbx.*GPglobal+ccx.*GPglobal.*GPglobal+ddx.*GPglobal.*GPglobal.*GPglobal;
    beta = aay+bby.*GPglobal+ccy.*GPglobal.*GPglobal+ddy.*GPglobal.*GPglobal.*GPglobal;
    deta = bbx+2*ccx.*GPglobal+3*ddx.*GPglobal.*GPglobal;
    dbeta = bby+2*ccy.*GPglobal+3*ddy.*GPglobal.*GPglobal;
    
    %chabge sign for my convention
    nx = dbeta./sqrt(deta.*deta+dbeta.*dbeta);
    ny = -deta./sqrt(deta.*deta+dbeta.*dbeta);
    
%     figure
%     plot(eta,nx)
%     hold on
%     plot(eta,ny,'r')
%     hold off
    
%     figure
%     plot(eta,beta,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off

%     h = zeros(1,el*6);
%     h0 = zeros(el,1);
%     h1 = zeros(el,1);

    %metrics term
%     for i = 1:el
%         h(1+6*(i-1):6*i) = sqrt(deta(1+6*(i-1):6*i).*deta(1+6*(i-1):6*i)+dbeta(1+6*(i-1):6*i).*dbeta(1+6*(i-1):6*i));
%         h0(i) = sqrt(bx(i)*bx(i)+by(i)*by(i));
%         h1(i) =  sqrt((bx(i)+2*cx(i)+3*dx(i))*(bx(i)+2*cx(i)+3*dx(i))+(by(i)+2*cy(i)+3*dy(i))*(by(i)+2*cy(i)+3*dy(i)));
%     end
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    h0 = sqrt(bx.*bx+by.*by);
    h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    
%     figure
%     plot(x(1:end-1),h0,'o-')
%     hold on
%     plot(x(2:end),h1,'r-o')
%     plot(eta,h,'k-o')
%     hold off

%     figure
%     plot(h,'o-')
    
    for i = 1:N
        for k = 1:6*el
            
            [SXX(i,k),SXY(i,k),SYX(i,k),SYY(i,k),QXXX(i,k),QXXY(i,k),QXYX(i,k),QXYY(i,k),QYXX(i,k),QYXY(i,k),QYYX(i,k),QYYY(i,k),PXX(i,k),PXY(i,k),PYX(i,k),PYY(i,k) ,Iaxis(i,k)] = sgf_ax_fs (2,eta(k),beta(k),X0(i),Y0(i));
            
            SXX(i,k) = SXX(i,k)*h(k);
            SXY(i,k) = SXY(i,k)*h(k);
            SYX(i,k) = SYX(i,k)*h(k);
            SYY(i,k) = SYY(i,k)*h(k);
            QXXX(i,k) = QXXX(i,k)*h(k);
            QXXY(i,k) = QXXY(i,k)*h(k);
            QXYX(i,k) = QXYX(i,k)*h(k);
            QXYY(i,k) = QXYY(i,k)*h(k);
            QYXX(i,k) = QYXX(i,k)*h(k);
            QYXY(i,k) = QYXY(i,k)*h(k);
            QYYX(i,k) = QYYX(i,k)*h(k);
            QYYY(i,k) = QYYY(i,k)*h(k);
            
        end
    end
    
%     figure
%     surf(SYY)
    
    %singularity treatment for linear element REPLACE WITH VECTORIAL
    %PRODUCT %%%%%
     for i = 2:N-2

%        SXX(i,1+6*(i-1):6*i) = SXX(i,1+6*(i-1):6*i) + 2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i)));
%        SXX(i+1,1+6*(i-1):6*i) = SXX(i+1,1+6*(i-1):6*i) + 2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i)));
              
       SYY(i,1+6*(i-1):6*i) = SYY(i,1+6*(i-1):6*i) + 2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)).*h(1+6*(i-1):6*i)...
           - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i)));
       SYY(i+1,1+6*(i-1):6*i) = SYY(i+1,1+6*(i-1):6*i) + 2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)).*h(1+6*(i-1):6*i)...
           - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i)));

    end
    
%        SXX(2,1:6) = SXX(2,1:6) + 2*log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)).*h(1:6)...
%            - 2*(log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1)));
%        SXX(N-1,1+6*(N-2):6*(N-1)) = SXX(N-1,1+6*(N-2):6*(N-1)) + 2*log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)).*h(1+6*(N-2):6*(N-1))...
%            - 2*(log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1)));
       
       %figure
       %surf(SXX)
    
       SYY(2,1:6) = SYY(2,1:6) + 2*log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)).*h(1:6)...
           - 2*(log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1)));
       SYY(N-1,1+6*(N-2):6*(N-1)) = SYY(N-1,1+6*(N-2):6*(N-1)) + 2*log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)).*h(1+6*(N-2):6*(N-1))...
           - 2*(log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1)));
       
%        %figure
%        %surf(SYY)
    
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
        for k = 1:el
            GXX(i,k) = GXX(i,k) + SXX(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW';
            GXX(i,k+1) = SXX(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW';
            
            GXY(i,k) = GXY(i,k) + SXY(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW';
            GXY(i,k+1) = SXY(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW';
            
            GYX(i,k) = GYX(i,k) + SYX(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW';
            GYX(i,k+1) = SYX(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW';
            
            GYY(i,k) = GYY(i,k) + SYY(i,1+6*(k-1):6*k).*phia(1+6*(k-1):6*k)*GW';
            GYY(i,k+1) = SYY(i,1+6*(k-1):6*k).*phib(1+6*(k-1):6*k)*GW';
            
            A11(i,k) = A11(i,k) + (QXXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QXXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phia(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            A11(i,k+1) = (QXXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QXXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6)).*phib(1+6*(k-1):6*k)...
                .*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            
            A12(i,k) = A12(i,k) + (QXYX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QXYY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phia(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            A12(i,k+1) = (QXYX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QXYY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phib(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            
            A21(i,k) = A21(i,k) + (QYXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QYXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phia(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            A21(i,k+1) = (QYXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QYXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phib(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            
            A22(i,k) = A22(i,k) + (QYYX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QYYY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phia(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            A22(i,k+1) = (QYYX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QYYY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*phib(1+6*(k-1):6*k).*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            
            T1(i,k) = (QXXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QXXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            T2(i,k) = (QYXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + QYXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            
            D1(i,k) = (PXX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + PXY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
            D2(i,k) = (PYX(i,1+6*(k-1):6*k).*nx(6*(k-1)+1:k*6) + PYY(i,1+6*(k-1):6*k).*ny(6*(k-1)+1:k*6))...
                .*sqrt(deta(1+6*(k-1):6*k).*deta(1+6*(k-1):6*k)+dbeta(1+6*(k-1):6*k).*dbeta(1+6*(k-1):6*k))*GW';
        end
    end
    
    %singularity treatment for linear curved element REPLACE WITH VECTORIAL
    %PRODUCT
    for i = 2:N-2
        
       GXX(i,i) = GXX(i,i) + 5/2*h0(i);
       GXX(i,i+1) = GXX(i,i+1) + 2*h0(i);
       GXX(i+1,i) = GXX(i+1,i) + 2*h1(i);
       GXX(i+1,i+1) = GXX(i+1,i+1) + 5/2*h1(i);
       
       GYY(i,i) = GYY(i,i) + 5/2*h0(i);
       GYY(i,i+1) = GYY(i,i+1) + 2*h0(i);
       GYY(i+1,i) = GYY(i+1,i) + 2*h1(i);
       GYY(i+1,i+1) = GYY(i+1,i+1) + 5/2*h1(i);
       
    end
    
%       figure
%       surf(GXX)
    
%       GXX(2,1) = GXX(2,1) + 2*h1(1);
%       GXX(2,2) = GXX(2,2) + 5/2*h1(1);
%       GXX(N-1,N-1) = GXX(N-1,N-1) + 5/2*h0(el);
%       GXX(N-1,N) = GXX(N-1,N) + 2*h0(el);
      
%       figure
%       surf(GXX)
    
      GYY(2,1) = GYY(2,1) + 2*h1(1);
      GYY(2,2) = GYY(2,2) + 5/2*h1(1);
      GYY(N-1,N-1) = GYY(N-1,N-1) + 5/2*h0(el);
      GYY(N-1,N) = GYY(N-1,N) + 2*h0(el);

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