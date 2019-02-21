%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2,Iaxis] = ...
    computeGT_only_2interface(x,y,drop1,ax,ay,bx,by,cx,cy,dx,dy)

    %tic

    walls = 0;

    %number of singularities and element
    N = numel(x);
    el = N-2;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    T1 = zeros(N,el);
    T2 = zeros(N,el);
    D1 = zeros(N,el);
    D2 = zeros(N,el);
    
    %gauss points and weigths
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
    GWt = GW/2; %careful!!!!!!!!!!!!!!!!!!!!
    
    %hat functions
    phia = 1-GPt;
    phib = GPt;
    
    % point where the variable will be stored
    X0 = x;
    Y0 = y;

    %points where I perform gauss integration
    GPglobal = repmat(GPt,1,el);
    
    %global coordinated retrived on the splines
    eta = aax+bbx.*GPglobal+ccx.*GPglobal.*GPglobal+ddx.*GPglobal.*GPglobal.*GPglobal;
    beta = aay+bby.*GPglobal+ccy.*GPglobal.*GPglobal+ddy.*GPglobal.*GPglobal.*GPglobal;
    deta = bbx+2*ccx.*GPglobal+3*ddx.*GPglobal.*GPglobal;
    dbeta = bby+2*ccy.*GPglobal+3*ddy.*GPglobal.*GPglobal;
    
    %change sign for my convention
    nx = dbeta./sqrt(deta.*deta+dbeta.*dbeta);
    ny = -deta./sqrt(deta.*deta+dbeta.*dbeta);
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    h0 = sqrt(bx.*bx+by.*by);
    h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    manyh = repmat(h,N,1)';
    
%     figure
%     plot([globalX_walls eta],[globalY_walls beta],'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);
        
    %compute Green's functions for stresses and velocities
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,PXX,PXY,PYX,PYY ,Iaxis] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    SXX = manyh.*SXX(6*walls+1:end,:);
    SXY = manyh.*SXY(6*walls+1:end,:);
    SYX = manyh.*SYX(6*walls+1:end,:);
    SYY = manyh.*SYY(6*walls+1:end,:);
    QXXX = manyh.*QXXX(6*walls+1:end,:);
    QXXY = manyh.*QXXY(6*walls+1:end,:);
    QXYX = manyh.*QXYX(6*walls+1:end,:);
    QXYY = manyh.*QXYY(6*walls+1:end,:);
    QYXX = manyh.*QYXX(6*walls+1:end,:);
    QYXY = manyh.*QYXY(6*walls+1:end,:);
    QYYX = manyh.*QYYX(6*walls+1:end,:);
    QYYY = manyh.*QYYY(6*walls+1:end,:);
    PXX = manyh.*PXX(6*walls+1:end,:);
    PXY = manyh.*PXY(6*walls+1:end,:);
    PYX = manyh.*PYX(6*walls+1:end,:);
    PYY = manyh.*PYY(6*walls+1:end,:);
    
    % singularity treatment spline elem
    %FIRST INTERFACE
    for i = 2:drop1-1

       SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-1-walls):6*(i-walls))-X0(i)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i)).^2)./GPt).*h(1+6*(i-walls-1):6*(i-walls)) + log(GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h0(i-walls))))';
       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+1)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-walls-1):6*(i-walls)) + log(1-GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h1(i-walls))))';
              
       SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-1-walls):6*(i-walls))-X0(i)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i)).^2)./GPt).*h(1+6*(i-walls-1):6*(i-walls)) + log(GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h0(i-walls))))';
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+1)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-walls-1):6*(i-walls)) + log(1-GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h1(i-walls))))';

    end
    
    %SECOND INTERFACE
    for i = drop1+2:el-1

       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i+1)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-1-walls):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i+1)).^2)./GPt).*h(1+6*(i-walls-1):6*(i-walls)) + log(GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h0(i-walls))))';
       SXX(1+6*(i-1):6*i,i+2) = SXX(1+6*(i-1):6*i,i+2) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+2)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+2)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+2)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+2)).^2)./(1-GPt)).*h(1+6*(i-walls-1):6*(i-walls)) + log(1-GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h1(i-walls))))';
              
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i+1)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-1-walls):6*(i-walls))-X0(i+1)).^2+(beta(1+6*(i-1-walls):6*(i-walls))-Y0(i+1)).^2)./GPt).*h(1+6*(i-walls-1):6*(i-walls)) + log(GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h0(i-walls))))';
       SYY(1+6*(i-1):6*i,i+2) = SYY(1+6*(i-1):6*i,i+2) + (2*log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+2)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+2)).^2)).*h(1+6*(i-walls-1):6*(i-walls))...
           - 2*(log(sqrt((eta(1+6*(i-walls-1):6*(i-walls))-X0(i+2)).^2+(beta(1+6*(i-walls-1):6*(i-walls))-Y0(i+2)).^2)./(1-GPt)).*h(1+6*(i-walls-1):6*(i-walls)) + log(1-GPt).*(h(1+6*(i-walls-1):6*(i-walls))-h1(i-walls))))';

    end
    
       %this because I don;t treat the point on the axis
       SXX(6*walls+1:6*walls+6,walls+2) = SXX(6*walls+1:6*walls+6,walls+2) + (2*log(sqrt((eta(1:6)-X0(walls+2)).^2+(beta(1:6)-Y0(walls+2)).^2)).*h(1:6)...
           - 2*(log(sqrt((eta(1:6)-X0(walls+2)).^2+(beta(1:6)-Y0(walls+2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
       SXX(6*(walls+drop1)+1:6*(walls+drop1)+6,walls+drop1+3) = SXX(6*(walls+drop1)+1:6*(walls+drop1)+6,walls+drop1+3) + (2*log(sqrt((eta(6*drop1+1:6*drop1+6)-X0(walls+drop1+3)).^2+(beta(6*drop1+1:6*drop1+6)-Y0(walls+drop1+3)).^2)).*h(6*drop1+1:6*drop1+6)...
          - 2*(log(sqrt((eta(6*drop1+1:6*drop1+6)-X0(walls+drop1+3)).^2+(beta(6*drop1+1:6*drop1+6)-Y0(walls+drop1+3)).^2)./(1-GPt)).*h(6*drop1+1:6*drop1+6) + log(1-GPt).*(h(6*drop1+1:6*drop1+6)-h1(drop1+1))))';
       SXX(1+6*(walls+drop1-1):6*(walls+drop1),walls+drop1) = SXX(1+6*(walls+drop1-1):6*(walls+drop1),walls+drop1) + (2*log(sqrt((eta(1+6*(drop1-1):6*drop1)-X0(walls+drop1)).^2+(beta(1+6*(drop1-1):6*drop1)-Y0(walls+drop1)).^2)).*h(1+6*(drop1-1):6*drop1)...
           - 2*(log(sqrt((eta(1+6*(drop1-1):6*drop1)-X0(walls+drop1)).^2+(beta(1+6*(drop1-1):6*drop1)-Y0(walls+drop1)).^2)./GPt).*h(1+6*(drop1-1):6*drop1) + log(GPt).*(h(1+6*(drop1-1):6*drop1)-h0(drop1))))';
       SXX(1+6*(el-1):end,N-1) = SXX(1+6*(el-1):6*el,N-1) + (2*log(sqrt((eta(1+6*(el-walls-1):6*(el-walls))-X0(N-1)).^2+(beta(1+6*(el-walls-1):6*(el-walls))-Y0(N-1)).^2)).*h(1+6*(el-walls-1):6*(el-walls))...
          - 2*(log(sqrt((eta(1+6*(el-walls-1):6*(el-walls))-X0(N-1)).^2+(beta(1+6*(el-walls-1):6*(el-walls))-Y0(N-1)).^2)./GPt).*h(1+6*(el-walls-1):6*(el-walls)) + log(GPt).*(h(1+6*(el-walls-1):6*(el-walls))-h0(el-walls))))';
       
       SYY(6*walls+1:6*walls+6,walls+2) = SYY(6*walls+1:6*walls+6,walls+2) + (2*log(sqrt((eta(1:6)-X0(walls+2)).^2+(beta(1:6)-Y0(walls+2)).^2)).*h(1:6)...
           - 2*(log(sqrt((eta(1:6)-X0(walls+2)).^2+(beta(1:6)-Y0(walls+2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
       SYY(6*(walls+drop1)+1:6*(walls+drop1)+6,walls+drop1+3) = SYY(6*(walls+drop1)+1:6*(walls+drop1)+6,walls+drop1+3) + (2*log(sqrt((eta(6*drop1+1:6*drop1+6)-X0(walls+drop1+3)).^2+(beta(6*drop1+1:6*drop1+6)-Y0(walls+drop1+3)).^2)).*h(6*drop1+1:6*drop1+6)...
          - 2*(log(sqrt((eta(6*drop1+1:6*drop1+6)-X0(walls+drop1+3)).^2+(beta(6*drop1+1:6*drop1+6)-Y0(walls+drop1+3)).^2)./(1-GPt)).*h(6*drop1+1:6*drop1+6) + log(1-GPt).*(h(6*drop1+1:6*drop1+6)-h1(drop1+1))))';
       SYY(1+6*(walls+drop1-1):6*(walls+drop1),walls+drop1) = SYY(1+6*(walls+drop1-1):6*(walls+drop1),walls+drop1) + (2*log(sqrt((eta(1+6*(drop1-1):6*drop1)-X0(walls+drop1)).^2+(beta(1+6*(drop1-1):6*drop1)-Y0(walls+drop1)).^2)).*h(1+6*(drop1-1):6*drop1)...
           - 2*(log(sqrt((eta(1+6*(drop1-1):6*drop1)-X0(walls+drop1)).^2+(beta(1+6*(drop1-1):6*drop1)-Y0(walls+drop1)).^2)./GPt).*h(1+6*(drop1-1):6*drop1) + log(GPt).*(h(1+6*(drop1-1):6*drop1)-h0(drop1))))';
       SYY(1+6*(el-1):end,N-1) = SYY(1+6*(el-1):6*el,N-1) + (2*log(sqrt((eta(1+6*(el-walls-1):6*(el-walls))-X0(N-1)).^2+(beta(1+6*(el-walls-1):6*(el-walls))-Y0(N-1)).^2)).*h(1+6*(el-walls-1):6*(el-walls))...
          - 2*(log(sqrt((eta(1+6*(el-walls-1):6*(el-walls))-X0(N-1)).^2+(beta(1+6*(el-walls-1):6*(el-walls))-Y0(N-1)).^2)./GPt).*h(1+6*(el-walls-1):6*(el-walls)) + log(GPt).*(h(1+6*(el-walls-1):6*(el-walls))-h0(el-walls))))';
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR SPLINE ELEMENTS INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %moltiplicate normals
    R1 = repmat(nx,N,1)';
    R2 = repmat(ny,N,1)';
    
    %for integration
    manyGW = repmat(GWt',el-walls,N);
    PHIA = repmat(phia',el-walls,N);
    PHIB = repmat(phib',el-walls,N);
            
    intSXXa = cumsum(SXX(1+6*walls:end,:).*manyGW.*PHIA);
    intSXXb = cumsum(SXX(1+6*walls:end,:).*manyGW.*PHIB);
    intSXYa = cumsum(SXY(1+6*walls:end,:).*manyGW.*PHIA);
    intSXYb = cumsum(SXY(1+6*walls:end,:).*manyGW.*PHIB);
    intSYXa = cumsum(SYX(1+6*walls:end,:).*manyGW.*PHIA);
    intSYXb = cumsum(SYX(1+6*walls:end,:).*manyGW.*PHIB);
    intSYYa = cumsum(SYY(1+6*walls:end,:).*manyGW.*PHIA);
    intSYYb = cumsum(SYY(1+6*walls:end,:).*manyGW.*PHIB);
    intQXXXa = cumsum(QXXX(1+6*walls:end,:).*manyGW.*PHIA.*R1);
    intQXXXb = cumsum(QXXX(1+6*walls:end,:).*manyGW.*PHIB.*R1);
    intQXXYa = cumsum(QXXY(1+6*walls:end,:).*manyGW.*PHIA.*R2);
    intQXXYb = cumsum(QXXY(1+6*walls:end,:).*manyGW.*PHIB.*R2);
    intQXYXa = cumsum(QXYX(1+6*walls:end,:).*manyGW.*PHIA.*R1);
    intQXYXb = cumsum(QXYX(1+6*walls:end,:).*manyGW.*PHIB.*R1);
    intQXYYa = cumsum(QXYY(1+6*walls:end,:).*manyGW.*PHIA.*R2);
    intQXYYb = cumsum(QXYY(1+6*walls:end,:).*manyGW.*PHIB.*R2);
    intQYXXa = cumsum(QYXX(1+6*walls:end,:).*manyGW.*PHIA.*R1);
    intQYXXb = cumsum(QYXX(1+6*walls:end,:).*manyGW.*PHIB.*R1);
    intQYXYa = cumsum(QYXY(1+6*walls:end,:).*manyGW.*PHIA.*R2);
    intQYXYb = cumsum(QYXY(1+6*walls:end,:).*manyGW.*PHIB.*R2);
    intQYYXa = cumsum(QYYX(1+6*walls:end,:).*manyGW.*PHIA.*R1);
    intQYYXb = cumsum(QYYX(1+6*walls:end,:).*manyGW.*PHIB.*R1);
    intQYYYa = cumsum(QYYY(1+6*walls:end,:).*manyGW.*PHIA.*R2);
    intQYYYb = cumsum(QYYY(1+6*walls:end,:).*manyGW.*PHIB.*R2);
    intQXXX_2nd = cumsum(QXXX(1+6*walls:end,:).*manyGW.*R1);
    intQXXY_2nd = cumsum(QXXY(1+6*walls:end,:).*manyGW.*R2);
    intQYXX_2nd = cumsum(QYXX(1+6*walls:end,:).*manyGW.*R1);
    intQYXY_2nd = cumsum(QYXY(1+6*walls:end,:).*manyGW.*R2);
    intPXX_2nd = cumsum(PXX(1+6*walls:end,:).*manyGW.*R1);
    intPXY_2nd = cumsum(PXY(1+6*walls:end,:).*manyGW.*R2);
    intPYX_2nd = cumsum(PYX(1+6*walls:end,:).*manyGW.*R1);
    intPYY_2nd = cumsum(PYY(1+6*walls:end,:).*manyGW.*R2);
    
    subSXXa(2:el-walls,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el-walls,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el-walls,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el-walls,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:el-walls,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:el-walls,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:el-walls,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el-walls,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el-walls,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el-walls,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el-walls,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el-walls,:) = intQXXYb(6:6:end-6,:);
    subQXYXa(2:el-walls,:) = intQXYXa(6:6:end-6,:);
    subQXYXb(2:el-walls,:) = intQXYXb(6:6:end-6,:);
    subQXYYa(2:el-walls,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el-walls,:) = intQXYYb(6:6:end-6,:);
    subQYXXa(2:el-walls,:) = intQYXXa(6:6:end-6,:);
    subQYXXb(2:el-walls,:) = intQYXXb(6:6:end-6,:);
    subQYXYa(2:el-walls,:) = intQYXYa(6:6:end-6,:);
    subQYXYb(2:el-walls,:) = intQYXYb(6:6:end-6,:);
    subQYYXa(2:el-walls,:) = intQYYXa(6:6:end-6,:);
    subQYYXb(2:el-walls,:) = intQYYXb(6:6:end-6,:);
    subQYYYa(2:el-walls,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el-walls,:) = intQYYYb(6:6:end-6,:);
    subQXXX_2nd(2:el-walls,:) = intQXXX_2nd(6:6:end-6,:);
    subQXXY_2nd(2:el-walls,:) = intQXXY_2nd(6:6:end-6,:);
    subQYXX_2nd(2:el-walls,:) = intQYXX_2nd(6:6:end-6,:);
    subQYXY_2nd(2:el-walls,:) = intQYXY_2nd(6:6:end-6,:);
    subPXX_2nd(2:el-walls,:) = intPXX_2nd(6:6:end-6,:);
    subPXY_2nd(2:el-walls,:) = intPXY_2nd(6:6:end-6,:);
    subPYX_2nd(2:el-walls,:) = intPYX_2nd(6:6:end-6,:);
    subPYY_2nd(2:el-walls,:) = intPYY_2nd(6:6:end-6,:);
    
    %compute integral
    GXX(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intSXXa(6:6:end,:)-subSXXa)';
    GXX(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        GXX(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intSXXb(6:6:end,:)-subSXXb)';
    GXY(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intSXYa(6:6:end,:)-subSXYa)';
    GXY(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        GXY(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intSXYb(6:6:end,:)-subSXYb)';
    GYX(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intSYXa(6:6:end,:)-subSYXa)';
    GYX(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        GYX(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intSYXb(6:6:end,:)-subSYXb)';
    GYY(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intSYYa(6:6:end,:)-subSYYa)';
    GYY(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        GYY(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intSYYb(6:6:end,:)-subSYYb)';
    
    A11(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intQXXXa(6:6:end,:)-subQXXXa + intQXXYa(6:6:end,:)-subQXXYa)';
    A11(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        A11(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intQXXXb(6:6:end,:)-subQXXXb + intQXXYb(6:6:end,:)-subQXXYb)';
    A12(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intQXYXa(6:6:end,:)-subQXYXa + intQXYYa(6:6:end,:)-subQXYYa)';
    A12(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        A12(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intQXYXb(6:6:end,:)-subQXYXb + intQXYYb(6:6:end,:)-subQXYYb)';
    A21(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intQYXXa(6:6:end,:)-subQYXXa + intQYXYa(6:6:end,:)-subQYXYa)';
    A21(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        A21(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intQYXXb(6:6:end,:)-subQYXXb + intQYXYb(6:6:end,:)-subQYXYb)';
    A22(:,[walls+1:walls+drop1 walls+drop1+2:N-1]) = (intQYYXa(6:6:end,:)-subQYYXa + intQYYYa(6:6:end,:)-subQYYYa)';
    A22(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) = ...
        A22(:,[walls+2:walls+drop1+1 walls+drop1+3:N]) + (intQYYXb(6:6:end,:)-subQYYXb + intQYYYb(6:6:end,:)-subQYYYb)';
    
    T1(:,walls+1:end) = (intQXXX_2nd(6:6:end,:)-subQXXX_2nd + intQXXY_2nd(6:6:end,:)-subQXXY_2nd)';
    T2(:,walls+1:end) = (intQYXX_2nd(6:6:end,:)-subQYXX_2nd + intQYXY_2nd(6:6:end,:)-subQYXY_2nd)';
    
    D1(:,walls+1:end) = (intPXX_2nd(6:6:end,:)-subPXX_2nd + intPXY_2nd(6:6:end,:)-subPXY_2nd)';
    D2(:,walls+1:end) = (intPYX_2nd(6:6:end,:)-subPYX_2nd + intPYY_2nd(6:6:end,:)-subPYY_2nd)';
    
    %singularity treatment spline element
    
    %FIRST INTERFACE
    GXX(walls+2:walls+drop1-1,walls+2:walls+drop1-1) =...
        GXX(walls+2:walls+drop1-1,walls+2:walls+drop1-1) +...
        diag(3/2*h0(2:drop1-1));
    GXX(walls+2:walls+drop1-1,walls+3:walls+drop1) =...
        GXX(walls+2:walls+drop1-1,walls+3:walls+drop1) +...
        diag(1/2*h0(2:drop1-1));
    GXX(walls+3:walls+drop1,walls+2:walls+drop1-1) =...
        GXX(walls+3:walls+drop1,walls+2:walls+drop1-1) +...
        diag(1/2*h1(2:drop1-1));
    GXX(walls+3:walls+drop1,walls+3:walls+drop1) =...
        GXX(walls+3:walls+drop1,walls+3:walls+drop1) +...
        diag(3/2*h1(2:drop1-1));
    
    %SECOND INTERFACE
    GXX(walls+drop1+3:N-2,walls+drop1+3:N-2) =...
        GXX(walls+drop1+3:N-2,walls+drop1+3:N-2) +...
        diag(3/2*h0(drop1+2:end-1));
    GXX(walls+drop1+3:N-2,walls+drop1+4:N-1) =...
        GXX(walls+drop1+3:N-2,walls+drop1+4:N-1) +...
        diag(1/2*h0(drop1+2:end-1));
    GXX(walls+drop1+4:N-1,walls+drop1+3:N-2) =...
        GXX(walls+drop1+4:N-1,walls+drop1+3:N-2) +...
        diag(1/2*h1(drop1+2:end-1));
    GXX(walls+drop1+4:N-1,walls+drop1+4:N-1) =...
        GXX(walls+drop1+4:N-1,walls+drop1+4:N-1) +...
        diag(3/2*h1(drop1+2:end-1));
    
    %FIRST INTERFACE
    GYY(walls+2:walls+drop1-1,walls+2:walls+drop1-1) =...
        GYY(walls+2:walls+drop1-1,walls+2:walls+drop1-1) +...
        diag(3/2*h0(2:drop1-1));
    GYY(walls+2:walls+drop1-1,walls+3:walls+drop1) =...
        GYY(walls+2:walls+drop1-1,walls+3:walls+drop1) +...
        diag(1/2*h0(2:drop1-1));
    GYY(walls+3:walls+drop1,walls+2:walls+drop1-1) =...
        GYY(walls+3:walls+drop1,walls+2:walls+drop1-1) +...
        diag(1/2*h1(2:drop1-1));
    GYY(walls+3:walls+drop1,walls+3:walls+drop1) =...
        GYY(walls+3:walls+drop1,walls+3:walls+drop1) +...
        diag(3/2*h1(2:drop1-1));
    
    %SECOND INTERFACE
    GYY(walls+drop1+3:N-2,walls+drop1+3:N-2) =...
        GYY(walls+drop1+3:N-2,walls+drop1+3:N-2) +...
        diag(3/2*h0(drop1+2:end-1));
    GYY(walls+drop1+3:N-2,walls+drop1+4:N-1) =...
        GYY(walls+drop1+3:N-2,walls+drop1+4:N-1) +...
        diag(1/2*h0(drop1+2:end-1));
    GYY(walls+drop1+4:N-1,walls+drop1+3:N-2) =...
        GYY(walls+drop1+4:N-1,walls+drop1+3:N-2) +...
        diag(1/2*h1(drop1+2:end-1));
    GYY(walls+drop1+4:N-1,walls+drop1+4:N-1) =...
        GYY(walls+drop1+4:N-1,walls+drop1+4:N-1) +...
        diag(3/2*h1(drop1+2:end-1));
    
    %because singularity on the axis don't has to be treated
    
    %FIRST INTERFACE
    GXX(walls+2,walls+1) = ...
        GXX(walls+2,walls+1) + diag(1/2*h1(1));
    GXX(walls+2,walls+2) =...
        GXX(walls+2,walls+2) + diag(3/2*h1(1));
    GXX(walls+drop1,walls+drop1) = ...
        GXX(walls+drop1,walls+drop1) + diag(3/2*h0(drop1));
    GXX(walls+drop1,walls+drop1+1) = ...
        GXX(walls+drop1,walls+drop1+1) + diag(1/2*h0(drop1));
    
    GYY(walls+2,walls+1) = ...
        GYY(walls+2,walls+1) + diag(1/2*h1(1));
    GYY(walls+2,walls+2) =...
        GYY(walls+2,walls+2) + diag(3/2*h1(1));
    GYY(walls+drop1,walls+drop1) = ...
        GYY(walls+drop1,walls+drop1) + diag(3/2*h0(drop1));
    GYY(walls+drop1,walls+drop1+1) = ...
        GYY(walls+drop1,walls+drop1+1) + diag(1/2*h0(drop1));
    
    %SECOND INTERFACE
    GXX(walls+drop1+3,walls+drop1+2) = ...
        GXX(walls+drop1+3,walls+drop1+2) + diag(1/2*h1(drop1+1));
    GXX(walls+drop1+3,walls+drop1+3) =...
        GXX(walls+drop1+3,walls+drop1+3) + diag(3/2*h1(drop1+1));
    GXX(N-1,N-1) = ...
        GXX(N-1,N-1) + diag(3/2*h0(el-walls));
    GXX(N-1,N) = ...
        GXX(N-1,N) + diag(1/2*h0(el-walls));
    
    GYY(walls+drop1+3,walls+drop1+2) = ...
        GYY(walls+drop1+3,walls+drop1+2) + diag(1/2*h1(drop1+1));
    GYY(walls+drop1+3,walls+drop1+3) =...
        GYY(walls+drop1+3,walls+drop1+3) + diag(3/2*h1(drop1+1));
    GYY(N-1,N-1) = ...
        GYY(N-1,N-1) + diag(3/2*h0(el-walls));
    GYY(N-1,N) = ...
        GYY(N-1,N) + diag(1/2*h0(el-walls));
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %T = toc

end