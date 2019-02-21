%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY,PXX,PXY,PYX,PYY ,Iaxis] = computeGT_arc(x,y)

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
    
    %compute the circ coeff of circles passing for free points (use symmetry for the ends points)
    [theta,dtheta,R,a,b,~] = CompArcElem(x,y);
    %[theta2,dtheta2,R2,a2,b2,~] = CompArcElemUp(x,y);
    %averaged values for the element
%     R = (R(1:end-1)+R(2:end))/2;
%     a = (a(1:end-1)+a(2:end))/2;
%     b = (b(1:end-1)+b(2:end))/2;
%     c = (c(1:end-1)+c(2:end))/2;

%     t1 = (x(2:end)-x(1:end-1))./sqrt((x(2:end)-x(1:end-1)).*x(2:end)-x(1:end-1)+(y(2:end)-y(1:end-1)).*(y(2:end)-y(1:end-1)));
%     t2 = (y(2:end)-y(1:end-1))./sqrt((x(2:end)-x(1:end-1)).*x(2:end)-x(1:end-1)+(y(2:end)-y(1:end-1)).*(y(2:end)-y(1:end-1)));
%     
%     scalar = [u' v'].*[t1' t2'];
%     ut = scalar(:,1)+scalar(:,2);
%     
%     theta = (ut>=0).*(theta1)'+(ut<0).*(theta2)';
%     dtheta = (ut>=0).*(dtheta1)'+(ut<0).*(dtheta2)';
%     R = (ut>=0).*(R1)'+(ut<0).*(R2)';
%     a = (ut>=0).*(a1)'+(ut<0).*(a2)';
%     b = (ut>=0).*(b1)'+(ut<0).*(b2)';
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored (0th order) PUT THEM ON THE ARC!!!
    X0 = R.*cos(theta)-a/2;
    Y0 = R.*sin(theta)-b/2;
    
%     %moltiplicate X0 and Y0
%     XX(1:6:6*N-5) = X0;
%     XX(2:6:6*N-4) = X0;
%     XX(3:6:6*N-3) = X0;
%     XX(4:6:6*N-2) = X0;
%     XX(5:6:6*N-1) = X0;
%     XX(6:6:6*N) = X0;
%     YY(1:6:6*N-5) = Y0;
%     YY(2:6:6*N-4) = Y0;
%     YY(3:6:6*N-3) = Y0;
%     YY(4:6:6*N-2) = Y0;
%     YY(5:6:6*N-1) = Y0;
%     YY(6:6:6*N) = Y0;
    
%     deltaX = abs(x(2:end)-x(1:end-1));
%     deltaY = abs(y(2:end)-y(1:end-1));
%     deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);

    TT(1:6:6*N-5) = theta;
    TT(2:6:6*N-4) = theta;
    TT(3:6:6*N-3) = theta;
    TT(4:6:6*N-2) = theta;
    TT(5:6:6*N-1) = theta;
    TT(6:6:6*N) = theta;
    RR(1:6:6*N-5) = R;
    RR(2:6:6*N-4) = R;
    RR(3:6:6*N-3) = R;
    RR(4:6:6*N-2) = R;
    RR(5:6:6*N-1) = R;
    RR(6:6:6*N) = R;
    aa(1:6:6*N-5) = a;
    aa(2:6:6*N-4) = a;
    aa(3:6:6*N-3) = a;
    aa(4:6:6*N-2) = a;
    aa(5:6:6*N-1) = a;
    aa(6:6:6*N) = a;
    bb(1:6:6*N-5) = b;
    bb(2:6:6*N-4) = b;
    bb(3:6:6*N-3) = b;
    bb(4:6:6*N-2) = b;
    bb(5:6:6*N-1) = b;
    bb(6:6:6*N) = b;
    
    
    % points where I perform gauss integration (projecting the gauss points on arcs)
    %every Gauss point
    GPtheta(1:6:6*N-5) = 0.5*GP(1)*dtheta;
    GPtheta(2:6:6*N-4) = 0.5*GP(2)*dtheta;
    GPtheta(3:6:6*N-3) = 0.5*GP(3)*dtheta;
    GPtheta(4:6:6*N-2) = 0.5*GP(4)*dtheta;
    GPtheta(5:6:6*N-1) = 0.5*GP(5)*dtheta;
    GPtheta(6:6:6*N) = 0.5*GP(6)*dtheta;
    globalTheta = TT+GPtheta;
    
    %compute normal and tangential vectors at every Gauss point for
    %integration
    %t = [-sin(globalTheta)' cos(globalTheta)'];
    n = [cos(globalTheta)' sin(globalTheta)'];
    
    
    globalX = RR.*cos(globalTheta)-aa/2;
    globalY = RR.*sin(globalTheta)-bb/2;
    
    arc_length = R.*dtheta;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
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
            GXX(i,k) = SXX(i,1+6*(k-1):6*k)*GW'*arc_length(k)/2;
            GXY(i,k) = SXY(i,1+6*(k-1):6*k)*GW'*arc_length(k)/2;
            GYX(i,k) = SYX(i,1+6*(k-1):6*k)*GW'*arc_length(k)/2;
            GYY(i,k) = SYY(i,1+6*(k-1):6*k)*GW'*arc_length(k)/2;
            TXXX(i,k) = QXXX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*arc_length(k)/2;
            TXXY(i,k) = QXXY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*arc_length(k)/2;
            TXYX(i,k) = QXYX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*arc_length(k)/2;
            TXYY(i,k) = QXYY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*arc_length(k)/2;
            TYXX(i,k) = QYXX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*arc_length(k)/2;
            TYXY(i,k) = QYXY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*arc_length(k)/2;
            TYYX(i,k) = QYYX(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,1))'*GW'*arc_length(k)/2;
            TYYY(i,k) = QYYY(i,1+6*(k-1):6*k).*(n(1+6*(k-1):6*k,2))'*GW'*arc_length(k)/2;
        end
        % singularity treatment
        %GXX(i,i) = GXX(i,i)+2*(-2*arc_length(i)/2*log(arc_length(i)/2)+3*arc_length(i)/2);
        %GYY(i,i) = GYY(i,i)+2*(-2*arc_length(i)/2*log(arc_length(i)/2)+3*arc_length(i)/2);
%         TXXX(i,i) = 0;
%         TXXY(i,i) = 0;
%         TXYX(i,i) = 0;
%         TXYY(i,i) = 0;
%         TYXX(i,i) = 0;
%         TYXY(i,i) = 0;
%         TYYX(i,i) = 0;
%         TYYY(i,i) = 0;
    end
    
    %T = toc

end