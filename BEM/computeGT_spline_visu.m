%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0,Y0)

    %tic

    %number of singularities
    N = numel(X0);
    %number of element
    el = numel(ax);
    
    GXX = zeros(N,el+1);
    GXY = zeros(N,el+1);
    GYX = zeros(N,el+1);
    GYY = zeros(N,el+1);
    A11 = zeros(N,el+1);
    A12 = zeros(N,el+1);
    A21 = zeros(N,el+1);
    A22 = zeros(N,el+1);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    aax = reshape((repmat(ax,6,1)),1,6*el);
    bbx = reshape((repmat(bx,6,1)),1,6*el);
    ccx = reshape((repmat(cx,6,1)),1,6*el);
    ddx = reshape((repmat(dx,6,1)),1,6*el);
    aay = reshape((repmat(ay,6,1)),1,6*el);
    bby = reshape((repmat(by,6,1)),1,6*el);
    ccy = reshape((repmat(cy,6,1)),1,6*el);
    ddy = reshape((repmat(dy,6,1)),1,6*el);
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;
    GW = GW/2;
    
    %hat functions
    phia = 1-GPt;
    phib = GPt;
    
    %points where I perform gauss integration
    GPglobal = repmat(GPt,1,el);
    
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
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    %h0 = sqrt(bx.*bx+by.*by);
    %h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    manyh = repmat(h,N,1)';
    
%     figure
%     plot(x(1:end-1),h0,'o-')
%     hold on
%     plot(x(2:end),h1,'r-o')
%     plot(eta,h,'k-o')
%     hold off

%     figure
%     plot(h,'o-')

    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);

    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,PXX,PXY,PYX,PYY] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    SXX = manyh.*SXX;
    SXY = manyh.*SXY;
    SYX = manyh.*SYX;
    SYY = manyh.*SYY;
    QXXX = manyh.*QXXX;
    QXXY = manyh.*QXXY;
    QXYX = manyh.*QXYX;
    QXYY = manyh.*QXYY;
    QYXX = manyh.*QYXX;
    QYXY = manyh.*QYXY;
    QYYX = manyh.*QYYX;
    QYYY = manyh.*QYYY;
    PXX = manyh.*PXX;
    PXY = manyh.*PXY;
    PYX = manyh.*PYX;
    PYY = manyh.*PYY;
    
%     figure
%     surf(SYY)
    
    %singularity treatment for linear element REPLACE WITH VECTORIAL
    %PRODUCT %%%%%
%      for i = 2:N-2
% 
%        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i))))';
%        SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i))))';
%               
%        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i))))';
%        SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i))))';
% 
%     end
%     
%        SXX(1:6,2) = SXX(1:6,2) + (2*log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)).*h(1:6)...
%            - 2*(log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
%        SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) + (2*log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)).*h(1+6*(N-2):6*(N-1))...
%            - 2*(log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1))))';
%     
%        SYY(1:6,2) = SYY(1:6,2) + (2*log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)).*h(1:6)...
%            - 2*(log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
%        SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) + (2*log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)).*h(1+6*(N-2):6*(N-1))...
%            - 2*(log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1))))';
       
    %moltiplicate normals
    R1 = repmat(nx,N,1)';
    R2 = repmat(ny,N,1)';
    
    %for integration
    manyGW = repmat(GW',el,N);
    PHIA = repmat(phia',el,N);
    PHIB = repmat(phib',el,N);
    
    intSXXa = cumsum(SXX.*manyGW.*PHIA);
    intSXXb = cumsum(SXX.*manyGW.*PHIB);
    intSXYa = cumsum(SXY.*manyGW.*PHIA);
    intSXYb = cumsum(SXY.*manyGW.*PHIB);
    intSYXa = cumsum(SYX.*manyGW.*PHIA);
    intSYXb = cumsum(SYX.*manyGW.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*PHIB);
    intQXXXa = cumsum(QXXX.*manyGW.*PHIA.*R1);
    intQXXXb = cumsum(QXXX.*manyGW.*PHIB.*R1);
    intQXXYa = cumsum(QXXY.*manyGW.*PHIA.*R2);
    intQXXYb = cumsum(QXXY.*manyGW.*PHIB.*R2);
    intQXYXa = cumsum(QXYX.*manyGW.*PHIA.*R1);
    intQXYXb = cumsum(QXYX.*manyGW.*PHIB.*R1);
    intQXYYa = cumsum(QXYY.*manyGW.*PHIA.*R2);
    intQXYYb = cumsum(QXYY.*manyGW.*PHIB.*R2);
    intQYXXa = cumsum(QYXX.*manyGW.*PHIA.*R1);
    intQYXXb = cumsum(QYXX.*manyGW.*PHIB.*R1);
    intQYXYa = cumsum(QYXY.*manyGW.*PHIA.*R2);
    intQYXYb = cumsum(QYXY.*manyGW.*PHIB.*R2);
    intQYYXa = cumsum(QYYX.*manyGW.*PHIA.*R1);
    intQYYXb = cumsum(QYYX.*manyGW.*PHIB.*R1);
    intQYYYa = cumsum(QYYY.*manyGW.*PHIA.*R2);
    intQYYYb = cumsum(QYYY.*manyGW.*PHIB.*R2);
    intQXXX = cumsum(QXXX.*manyGW.*R1);
    intQXXY = cumsum(QXXY.*manyGW.*R2);
    intQYXX = cumsum(QYXX.*manyGW.*R1);
    intQYXY = cumsum(QYXY.*manyGW.*R2);
    intPXX = cumsum(PXX.*manyGW.*R1);
    intPXY = cumsum(PXY.*manyGW.*R2);
    intPYX = cumsum(PYX.*manyGW.*R1);
    intPYY = cumsum(PYY.*manyGW.*R2);
    
    subSXXa(2:el,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:el,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:el,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:el,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el,:) = intQXXYb(6:6:end-6,:);
    subQXYXa(2:el,:) = intQXYXa(6:6:end-6,:);
    subQXYXb(2:el,:) = intQXYXb(6:6:end-6,:);
    subQXYYa(2:el,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el,:) = intQXYYb(6:6:end-6,:);
    subQYXXa(2:el,:) = intQYXXa(6:6:end-6,:);
    subQYXXb(2:el,:) = intQYXXb(6:6:end-6,:);
    subQYXYa(2:el,:) = intQYXYa(6:6:end-6,:);
    subQYXYb(2:el,:) = intQYXYb(6:6:end-6,:);
    subQYYXa(2:el,:) = intQYYXa(6:6:end-6,:);
    subQYYXb(2:el,:) = intQYYXb(6:6:end-6,:);
    subQYYYa(2:el,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el,:) = intQYYYb(6:6:end-6,:);
    subQXXX(2:el,:) = intQXXX(6:6:end-6,:);
    subQXXY(2:el,:) = intQXXY(6:6:end-6,:);
    subQYXX(2:el,:) = intQYXX(6:6:end-6,:);
    subQYXY(2:el,:) = intQYXY(6:6:end-6,:);
    subPXX(2:el,:) = intPXX(6:6:end-6,:);
    subPXY(2:el,:) = intPXY(6:6:end-6,:);
    subPYX(2:el,:) = intPYX(6:6:end-6,:);
    subPYY(2:el,:) = intPYY(6:6:end-6,:);
    
    %compute integral
    GXX(:,1:el) = (intSXXa(6:6:end,:)-subSXXa)';
    GXX(:,2:el+1) = GXX(:,2:el+1) + (intSXXb(6:6:end,:)-subSXXb)';
    GXY(:,1:el) = (intSXYa(6:6:end,:)-subSXYa)';
    GXY(:,2:el+1) = GXY(:,2:el+1) + (intSXYb(6:6:end,:)-subSXYb)';
    GYX(:,1:el) = (intSYXa(6:6:end,:)-subSYXa)';
    GYX(:,2:el+1) = GYX(:,2:el+1) + (intSYXb(6:6:end,:)-subSYXb)';
    GYY(:,1:el) = (intSYYa(6:6:end,:)-subSYYa)';
    GYY(:,2:el+1) = GYY(:,2:el+1) + (intSYYb(6:6:end,:)-subSYYb)';
    
    A11(:,1:el) = (intQXXXa(6:6:end,:)-subQXXXa + intQXXYa(6:6:end,:)-subQXXYa)';
    A11(:,2:el+1) = A11(:,2:el+1) + (intQXXXb(6:6:end,:)-subQXXXb + intQXXYb(6:6:end,:)-subQXXYb)';
    A12(:,1:el) = (intQXYXa(6:6:end,:)-subQXYXa + intQXYYa(6:6:end,:)-subQXYYa)';
    A12(:,2:el+1) = A12(:,2:el+1) + (intQXYXb(6:6:end,:)-subQXYXb + intQXYYb(6:6:end,:)-subQXYYb)';
    A21(:,1:el) = (intQYXXa(6:6:end,:)-subQYXXa + intQYXYa(6:6:end,:)-subQYXYa)';
    A21(:,2:el+1) = A21(:,2:el+1) + (intQYXXb(6:6:end,:)-subQYXXb + intQYXYb(6:6:end,:)-subQYXYb)';
    A22(:,1:el) = (intQYYXa(6:6:end,:)-subQYYXa + intQYYYa(6:6:end,:)-subQYYYa)';
    A22(:,2:el+1) = A22(:,2:el+1) + (intQYYXb(6:6:end,:)-subQYYXb + intQYYYb(6:6:end,:)-subQYYYb)';
    
    T1 = (intQXXX(6:6:end,:)-subQXXX + intQXXY(6:6:end,:)-subQXXY)';
    T2 = (intQYXX(6:6:end,:)-subQYXX + intQYXY(6:6:end,:)-subQYXY)';
    
    D1 = (intPXX(6:6:end,:)-subPXX + intPXY(6:6:end,:)-subPXY)';
    D2 = (intPYX(6:6:end,:)-subPYX + intPYY(6:6:end,:)-subPYY)';
    
    %singularity treatment
%     GXX(2:N-2,2:N-2) = GXX(2:N-2,2:N-2) + diag(3/2*h0(2:N-2));
%     GXX(2:N-2,3:N-1) = GXX(2:N-2,3:N-1) + diag(1/2*h0(2:N-2));
%     GXX(3:N-1,2:N-2) = GXX(3:N-1,2:N-2) + diag(1/2*h1(2:N-2));
%     GXX(3:N-1,3:N-1) = GXX(3:N-1,3:N-1) + diag(3/2*h1(2:N-2));
%     
%     GYY(2:N-2,2:N-2) = GYY(2:N-2,2:N-2) + diag(3/2*h0(2:N-2));
%     GYY(2:N-2,3:N-1) = GYY(2:N-2,3:N-1) + diag(1/2*h0(2:N-2));
%     GYY(3:N-1,2:N-2) = GYY(3:N-1,2:N-2) + diag(1/2*h1(2:N-2));
%     GYY(3:N-1,3:N-1) = GYY(3:N-1,3:N-1) + diag(3/2*h1(2:N-2));
%     
%     %because singularity on the axis don't has to be treated
%     GXX(2,1) = GXX(2,1) + 1/2*h1(1);
%     GXX(2,2) = GXX(2,2) + 3/2*h1(1);
%     GXX(N-1,N-1) = GXX(N-1,N-1) + 3/2*h0(el);
%     GXX(N-1,N) = GXX(N-1,N) + 1/2*h0(el);
%     
%     GYY(2,1) = GYY(2,1) + 1/2*h1(1);
%     GYY(2,2) = GYY(2,2) + 3/2*h1(1);
%     GYY(N-1,N-1) = GYY(N-1,N-1) + 3/2*h0(el);
%     GYY(N-1,N) = GYY(N-1,N) + 1/2*h0(el);
    
    %T = toc

end