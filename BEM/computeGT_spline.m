%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY,PXX,PXY,PYX,PYY ,Iaxis] = computeGT_spline(x,y)

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
    
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = my_spline (x, y);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored (0th order) PUT THEM ON THE ARC!!!
    X0 = ax+bx*0.5+cx*0.5^2+dx*0.5^3;
    Y0 = ay+by*0.5+cy*0.5^2+dy*0.5^3;
    
    aax(1:6:6*N-5) = ax;
    aax(2:6:6*N-4) = ax;
    aax(3:6:6*N-3) = ax;
    aax(4:6:6*N-2) = ax;
    aax(5:6:6*N-1) = ax;
    aax(6:6:6*N) = ax;
    bbx(1:6:6*N-5) = bx;
    bbx(2:6:6*N-4) = bx;
    bbx(3:6:6*N-3) = bx;
    bbx(4:6:6*N-2) = bx;
    bbx(5:6:6*N-1) = bx;
    bbx(6:6:6*N) = bx;
    ccx(1:6:6*N-5) = cx;
    ccx(2:6:6*N-4) = cx;
    ccx(3:6:6*N-3) = cx;
    ccx(4:6:6*N-2) = cx;
    ccx(5:6:6*N-1) = cx;
    ccx(6:6:6*N) = cx;
    ddx(1:6:6*N-5) = dx;
    ddx(2:6:6*N-4) = dx;
    ddx(3:6:6*N-3) = dx;
    ddx(4:6:6*N-2) = dx;
    ddx(5:6:6*N-1) = dx;
    ddx(6:6:6*N) = dx;
    aay(1:6:6*N-5) = ay;
    aay(2:6:6*N-4) = ay;
    aay(3:6:6*N-3) = ay;
    aay(4:6:6*N-2) = ay;
    aay(5:6:6*N-1) = ay;
    aay(6:6:6*N) = ay;
    bby(1:6:6*N-5) = by;
    bby(2:6:6*N-4) = by;
    bby(3:6:6*N-3) = by;
    bby(4:6:6*N-2) = by;
    bby(5:6:6*N-1) = by;
    bby(6:6:6*N) = by;
    ccy(1:6:6*N-5) = cy;
    ccy(2:6:6*N-4) = cy;
    ccy(3:6:6*N-3) = cy;
    ccy(4:6:6*N-2) = cy;
    ccy(5:6:6*N-1) = cy;
    ccy(6:6:6*N) = cy;
    ddy(1:6:6*N-5) = dy;
    ddy(2:6:6*N-4) = dy;
    ddy(3:6:6*N-3) = dy;
    ddy(4:6:6*N-2) = dy;
    ddy(5:6:6*N-1) = dy;
    ddy(6:6:6*N) = dy;
    
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;    
    
    % points where I perform gauss integration
    GPglobal(1:6:6*N-5) = GPt(1);
    GPglobal(2:6:6*N-4) = GPt(2);
    GPglobal(3:6:6*N-3) = GPt(3);
    GPglobal(4:6:6*N-2) = GPt(4);
    GPglobal(5:6:6*N-1) = GPt(5);
    GPglobal(6:6:6*N) = GPt(6);
    
    globalX = aax+bbx.*GPglobal+ccx.*(GPglobal.*GPglobal)+ddx.*GPglobal.*(GPglobal.*GPglobal.*GPglobal);
    globalY = aay+bby.*GPglobal+ccy.*(GPglobal.*GPglobal)+ddy.*GPglobal.*(GPglobal.*GPglobal.*GPglobal);
    
    %compute normal and vectors at every Gauss point for
    %integration    
%     [~,nx,ny ] = normal_cubic_spline(x,y);
%     
%     %chabge sign for my convention
%     nx = -nx;
%     ny = -ny;
    
    %arc_length = R.*dtheta;
    el = zeros(numel(x)-1,1);
    for i = 1:numel(x)-1
        f = @(t) sqrt((bx(i)+2*cx(i)*t+3*dx(i)*t).^2+(by(i)+2*cy(i)*t+3*dy(i)*t).^2);
        el(i) = integral(f,0,1);
    end
    
    figure
    plot(globalX,globalY,'o-')
    hold on
    axis equal
    plot(X0,Y0,'og-')
    hold off
    
    for i = 1:N
        for k = 1:6*N
            [SXX(i,k),SXY(i,k),SYX(i,k),SYY(i,k),QXXX(i,k),QXXY(i,k),QXYX(i,k),QXYY(i,k),QYXX(i,k),QYXY(i,k),QYYX(i,k),QYYY(i,k),PXX(i,k),PXY(i,k),PYX(i,k),PYY(i,k) ,Iaxis(i,k)] = sgf_ax_fs (2,globalX(k),globalY(k),X0(i),Y0(i));
            
        end
        %singularity treatment CORREGGI PER CURVED ELEMENT
        %SXX(i,1+6*(i-1):6*i) = SXX(i,1+6*(i-1):6*i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))-1;
        %SYY(i,1+6*(i-1):6*i) = SYY(i,1+6*(i-1):6*i)+2*log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))-1;
    end
    
    
    for i = 1:N
        for k = 1:N
            GXX(i,k) = SXX(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            GXY(i,k) = SXY(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            GYX(i,k) = SYX(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            GYY(i,k) = SYY(i,1+6*(k-1):6*k)*GW'*el(k)/2;
%             TXXX(i,k) = QXXX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*el(k)/2;
%             TXXY(i,k) = QXXY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*el(k)/2;
%             TXYX(i,k) = QXYX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*el(k)/2;
%             TXYY(i,k) = QXYY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*el(k)/2;
%             TYXX(i,k) = QYXX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*el(k)/2;
%             TYXY(i,k) = QYXY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*el(k)/2;
%             TYYX(i,k) = QYYX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*el(k)/2;
%             TYYY(i,k) = QYYY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*el(k)/2;
            TXXX(i,k) = QXXX(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TXXY(i,k) = QXXY(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TXYX(i,k) = QXYX(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TXYY(i,k) = QXYY(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TYXX(i,k) = QYXX(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TYXY(i,k) = QYXY(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TYYX(i,k) = QYYX(i,1+6*(k-1):6*k)*GW'*el(k)/2;
            TYYY(i,k) = QYYY(i,1+6*(k-1):6*k)*GW'*el(k)/2;
        end
        % singularity treatment
        %GXX(i,i) = GXX(i,i)+2*(-2*arc_length(i)/2*log(arc_length(i)/2)+3*arc_length(i)/2);
        %GYY(i,i) = GYY(i,i)+2*(-2*arc_length(i)/2*log(arc_length(i)/2)+3*arc_length(i)/2);
    end
    
    %T = toc

end