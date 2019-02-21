%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,PI1,PI2] = computeKernelPressureStokesAxisConstSplineNoSing(x0,y0,ax,ay,bx,by,cx,cy,dx,dy)

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
    [LX,LY,NXX,NXY,NYX,NYY] = kernelPressureAxis(globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %mutuply kernnel of doublelayer times the normal
    LX = manyh.*LX;
    LY = manyh.*LY;
    N1 = manyh.*NXX.*nnx + manyh.*NXY.*nny;
    N2 = manyh.*NYX.*nnx + manyh.*NYY.*nny;
    
    %INTEGRATION
    manyGW = repmat(GW',elem,N);
    
    %INTEGRATION
    intLX = cumsum(LX.*manyGW);
    intLY = cumsum(LY.*manyGW);
    intN1 = cumsum(N1.*manyGW);
    intN2 = cumsum(N2.*manyGW);
    
    subLX(2:elem,:) = intLX(6:6:end-6,:);
    subLY(2:elem,:) = intLY(6:6:end-6,:);
    subN1(2:elem,:) = intN1(6:6:end-6,:);
    subN2(2:elem,:) = intN2(6:6:end-6,:);
    
    %compute integral single layer
    PX = (intLX(6:6:end,:)-subLX)';
    PY = (intLY(6:6:end,:)-subLY)';
    
    %compute integral double layer
    PI1 = (intN1(6:6:end,:)-subN1)';
    PI2 = (intN2(6:6:end,:)-subN2)';

end