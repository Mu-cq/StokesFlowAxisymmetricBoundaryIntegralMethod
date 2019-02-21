%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelStokesAxisConstSplineNoSing(x0,y0,ax,ay,bx,by,cx,cy,dx,dy)

    %number of singularities
    N = numel(x0);
    elem = numel(ax);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    aax = reshape((repmat(ax,6,1)),1,6*elem);
    bbx = reshape((repmat(bx,6,1)),1,6*elem);
    ccx = reshape((repmat(cx,6,1)),1,6*elem);
    ddx = reshape((repmat(dx,6,1)),1,6*elem);
    aay = reshape((repmat(ay,6,1)),1,6*elem);
    bby = reshape((repmat(by,6,1)),1,6*elem);
    ccy = reshape((repmat(cy,6,1)),1,6*elem);
    ddy = reshape((repmat(dy,6,1)),1,6*elem);    
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;
    GW = GW/2;

    % point where the variable will be stored
    X0 = x0';
    Y0 = y0';
    
    %points where I perform gauss integration
    GPglobal = repmat(GPt,1,elem);
    
    %global coordianted retrived on the splines
    eta = aax+bbx.*GPglobal+ccx.*GPglobal.*GPglobal+ddx.*GPglobal.*GPglobal.*GPglobal;
    beta = aay+bby.*GPglobal+ccy.*GPglobal.*GPglobal+ddy.*GPglobal.*GPglobal.*GPglobal;
    deta = bbx+2*ccx.*GPglobal+3*ddx.*GPglobal.*GPglobal;
    dbeta = bby+2*ccy.*GPglobal+3*ddy.*GPglobal.*GPglobal;
    
    %chabge sign for my convention
    nx = dbeta./sqrt(deta.*deta+dbeta.*dbeta);
    ny = -deta./sqrt(deta.*deta+dbeta.*dbeta);
    
    %moltiplicate normals
    nnx = repmat(nx,N,1)';
    nny = repmat(ny,N,1)';
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    manyh = repmat(h,N,1)';

    X0matr = repmat(X0,6*elem,1);
    Y0matr = repmat(Y0,6*elem,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);
        
    %compute gree function and the gradient
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
          sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %mutuply kernnel of doublelayer times the normal
    SXX = manyh.*SXX;
    SXY = manyh.*SXY;
    SYX = manyh.*SYX;
    SYY = manyh.*SYY;
    Q11 = manyh.*QXXX.*nnx + manyh.*QXXY.*nny;
    Q12 = manyh.*QXYX.*nnx + manyh.*QXYY.*nny;
    Q21 = manyh.*QYXX.*nnx + manyh.*QYXY.*nny;
    Q22 = manyh.*QYYX.*nnx + manyh.*QYYY.*nny;
    
    %INTEGRATION
    manyGW = repmat(GW',elem,N);
    
    %INTEGRATION
    intSXX = cumsum(SXX.*manyGW);
    intSXY = cumsum(SXY.*manyGW);
    intSYX = cumsum(SYX.*manyGW);
    intSYY = cumsum(SYY.*manyGW);
    intQ11 = cumsum(Q11.*manyGW);
    intQ12 = cumsum(Q12.*manyGW);
    intQ21 = cumsum(Q21.*manyGW);
    intQ22 = cumsum(Q22.*manyGW);
    
    subSXX(2:elem,:) = intSXX(6:6:end-6,:);
    subSXY(2:elem,:) = intSXY(6:6:end-6,:);
    subSYX(2:elem,:) = intSYX(6:6:end-6,:);
    subSYY(2:elem,:) = intSYY(6:6:end-6,:);
    subQ11(2:elem,:) = intQ11(6:6:end-6,:);
    subQ12(2:elem,:) = intQ12(6:6:end-6,:);
    subQ21(2:elem,:) = intQ21(6:6:end-6,:);
    subQ22(2:elem,:) = intQ22(6:6:end-6,:);
    
    %compute integral single layer
    GXX = (intSXX(6:6:end,:)-subSXX)';
    GXY = (intSXY(6:6:end,:)-subSXY)';
    GYX = (intSYX(6:6:end,:)-subSYX)';
    GYY = (intSYY(6:6:end,:)-subSYY)';
    
    %compute integral double layer
    A11 = (intQ11(6:6:end,:)-subQ11)';
    A12 = (intQ12(6:6:end,:)-subQ12)';
    A21 = (intQ21(6:6:end,:)-subQ21)';
    A22 = (intQ22(6:6:end,:)-subQ22)';

end