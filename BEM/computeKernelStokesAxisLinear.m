%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelStokesAxisLinear(x,y,x0,y0,indIn,indOut,ST,kernelFreeSpace,wallPos)

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
    %initialize
    GXX = zeros(N,elem+1);
    GXY = zeros(N,elem+1);
    GYX = zeros(N,elem+1);
    GYY = zeros(N,elem+1);
    A11 = zeros(N,elem+1);
    A12 = zeros(N,elem+1);
    A21 = zeros(N,elem+1);
    A22 = zeros(N,elem+1);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x0';
    Y0 = y0';
    
    %compute normal vectors
    no = sqrt(diff(x).^2+diff(y).^2);
    nx = diff(y)./no;
    ny = -diff(x)./no;
    nnx = repmat(nx,6,1);
    nny = repmat(ny,6,1);
    nnx = repmat(nnx(:),1,N);
    nny = repmat(nny(:),1,N);
    
    %moltiplicate first point of the elment
    tempx = repmat(x(1:end-1),6,1);
    XX = reshape(tempx,1,6*elem);
    tempy = repmat(y(1:end-1),6,1);
    YY = reshape(tempy,1,6*elem);
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %modify GP
    GP = GP+1;
    
    %every Gauss point
    GPX = repmat(GP,1,elem).*reshape((repmat(deltaX/2,6,1)),1,6*elem);
    GPY = repmat(GP,1,elem).*reshape((repmat(deltaY/2,6,1)),1,6*elem);
    
    phia = 1-GP/2;
    phib = GP/2;
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
    X0matr = repmat(X0,6*elem,1);
    Y0matr = repmat(Y0,6*elem,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %compute gree function and the gradient
    if kernelFreeSpace==1     %free space
          [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
            sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
    elseif kernelFreeSpace==2     %free space
          [SXX,SXY,SYX,SYY ...
            ,QXXX,QXXY,QXYX,QXYY ...
            ,QYXX,QYXY,QYYX,QYYY] = sgf_ax_w_vect (globalXmatr,globalYmatr,X0matr,Y0matr,wallPos);
    else
        error('Not implemented')
    end
    
    %mutuply kernnel of doublelayer times the normal
    Q11 = QXXX.*nnx + QXXY.*nny;
    Q12 = QXYX.*nnx + QXYY.*nny;
    Q21 = QYXX.*nnx + QYXY.*nny;
    Q22 = QYYX.*nnx + QYYY.*nny;
      
    %singular treatment
    if ST==1
        %error('not implmented')
        range = indIn:indOut;
        [SXX(:,range),SYY(:,range),anA,anB] = stokesSTlinear(SXX(:,range),SYY(:,range),deltaL,globalXmatr(:,range),globalYmatr(:,range),x0(range),y0(range));
    end
    
    %INTEGRATION
    temp = reshape(repmat(deltaL/2,6,1),1,6*elem);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',elem,N);
    PHIA = repmat(phia',elem,N);
    PHIB = repmat(phib',elem,N);
    
    %INTEGRATION
    intSXXa = cumsum(SXX.*manyGW.*manyDelta.*PHIA);
    intSXXb = cumsum(SXX.*manyGW.*manyDelta.*PHIB);
    intSXYa = cumsum(SXY.*manyGW.*manyDelta.*PHIA);
    intSXYb = cumsum(SXY.*manyGW.*manyDelta.*PHIB);
    intSYXa = cumsum(SYX.*manyGW.*manyDelta.*PHIA);
    intSYXb = cumsum(SYX.*manyGW.*manyDelta.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*manyDelta.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*manyDelta.*PHIB);
    intQ11a = cumsum(Q11.*manyGW.*manyDelta.*PHIA);
    intQ11b = cumsum(Q11.*manyGW.*manyDelta.*PHIB);
    intQ12a = cumsum(Q12.*manyGW.*manyDelta.*PHIA);
    intQ12b = cumsum(Q12.*manyGW.*manyDelta.*PHIB);
    intQ21a = cumsum(Q21.*manyGW.*manyDelta.*PHIA);
    intQ21b = cumsum(Q21.*manyGW.*manyDelta.*PHIB);
    intQ22a = cumsum(Q22.*manyGW.*manyDelta.*PHIA);
    intQ22b = cumsum(Q22.*manyGW.*manyDelta.*PHIB);
    
    subSXXa(2:elem,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:elem,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:elem,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:elem,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:elem,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:elem,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:elem,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:elem,:) = intSYYb(6:6:end-6,:);
    subQ11a(2:elem,:) = intQ11a(6:6:end-6,:);
    subQ11b(2:elem,:) = intQ11b(6:6:end-6,:);
    subQ12a(2:elem,:) = intQ12a(6:6:end-6,:);
    subQ12b(2:elem,:) = intQ12b(6:6:end-6,:);
    subQ21a(2:elem,:) = intQ21a(6:6:end-6,:);
    subQ21b(2:elem,:) = intQ21b(6:6:end-6,:);
    subQ22a(2:elem,:) = intQ22a(6:6:end-6,:);
    subQ22b(2:elem,:) = intQ22b(6:6:end-6,:);
    
    %compute integral
    GXX(:,1:elem) = (intSXXa(6:6:end,:)-subSXXa)'; %+ vertcat(zeros(1,N-1),diag(-deltaL(1:N-1).*log(deltaL(1:N-1))+deltaL(1:N-1)));
    GXX(:,2:elem+1) = GXX(:,2:elem+1) + intSXXb(6:6:end,:)'-subSXXb';
    GXY(:,1:elem) = intSXYa(6:6:end,:)'-subSXYa';
    GXY(:,2:elem+1) = GXY(:,2:elem+1) + intSXYb(6:6:end,:)'-subSXYb';
    GYX(:,1:elem) = intSYXa(6:6:end,:)'-subSYXa';
    GYX(:,2:elem+1) = GYX(:,2:elem+1) + intSYXb(6:6:end,:)'-subSYXb';
    GYY(:,1:elem) = intSYYa(6:6:end,:)'-subSYYa';
    GYY(:,2:elem+1) = GYY(:,2:elem+1) + intSYYb(6:6:end,:)'-subSYYb';
    
    A11(:,1:elem) = (intQ11a(6:6:end,:)-subQ11a)';
    A11(:,2:elem+1) = A11(:,2:elem+1) + (intQ11b(6:6:end,:)-subQ11b)';
    A12(:,1:elem) = (intQ12a(6:6:end,:)-subQ12a)';
    A12(:,2:elem+1) = A12(:,2:elem+1) + (intQ12b(6:6:end,:)-subQ12b)';
    A21(:,1:elem) = (intQ21a(6:6:end,:)-subQ21a)';
    A21(:,2:elem+1) = A21(:,2:elem+1) + (intQ21b(6:6:end,:)-subQ21b)';
    A22(:,1:elem) = (intQ22a(6:6:end,:)-subQ22a)';
    A22(:,2:elem+1) = A22(:,2:elem+1) + (intQ22b(6:6:end,:)-subQ22b)';
    
    %add anaytical integration of the singularity
    if ST==1
        
        N = numel(deltaL)+1;
        
        range1 = (2:N-2)+indIn(1)-1;
        range2 = (2:N-2)+indIn(1);
        
        %add analtyical integration
        GXX(range1,2:N-2) = GXX(range1,2:N-2) + anA;
        GXX(range1,3:N-1) = GXX(range1,3:N-1) + anB;
        GXX(range2,2:N-2) = GXX(range2,2:N-2) + anB;
        GXX(range2,3:N-1) = GXX(range2,3:N-1) + anA;
        
        GYY(range1,2:N-2) = GYY(range1,2:N-2) + anA;
        GYY(range1,3:N-1) = GYY(range1,3:N-1) + anB;
        GYY(range2,2:N-2) = GYY(range2,2:N-2) + anB;
        GYY(range2,3:N-1) = GYY(range2,3:N-1) + anA;
        
        %because singularity on the axis don't has to be treated
        GXX(range(2),1) = GXX(range(2),1) - deltaL(1)*log(deltaL(1))+deltaL(1);
        GXX(range(2),2) = GXX(range(2),2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
        GXX(range(end-1),end-1) = GXX(range(end-1),end-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
        GXX(range(end-1),end) = GXX(range(end-1),end) - deltaL(N-1)*log(deltaL(range(end-1)))+deltaL(N-1);
        
        GYY(range(2),1) = GYY(range(2),1) - deltaL(1)*log(deltaL(1))+deltaL(1);
        GYY(range(2),2) = GYY(range(2),2) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
        GYY(range(end-1),end-1) = GYY(range(end-1),end-1) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
        GYY(range(end-1),end) = GYY(range(end-1),end) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
    
    end

end