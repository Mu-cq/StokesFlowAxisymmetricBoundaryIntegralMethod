%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,PI1,PI2] = computeKernelPressureStokesAxisLinearSpline(x0,y0,ax,ay,bx,by,cx,cy,dx,dy)

    %number of singularities
    N = numel(x0);
    elem = numel(ax);
    
    %initialize
    PX = zeros(N,elem+1);
    PY = zeros(N,elem+1);
    PI1 = zeros(N,elem+1);
    PI2 = zeros(N,elem+1);
    
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
    
    phia = 1-GPt;
    phib = GPt;

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
    PHIA = repmat(phia',elem,N);
    PHIB = repmat(phib',elem,N);
    
    %INTEGRATION
    intLXa = cumsum(LX.*manyGW.*PHIA);
    intLXb = cumsum(LX.*manyGW.*PHIB);
    intLYa = cumsum(LY.*manyGW.*PHIA);
    intLYb = cumsum(LY.*manyGW.*PHIB);
    intN1a = cumsum(N1.*manyGW.*PHIA);
    intN1b = cumsum(N1.*manyGW.*PHIB);
    intN2a = cumsum(N2.*manyGW.*PHIA);
    intN2b = cumsum(N2.*manyGW.*PHIB);
    
    subLXa(2:elem,:) = intLXa(6:6:end-6,:);
    subLXb(2:elem,:) = intLXb(6:6:end-6,:);
    subLYa(2:elem,:) = intLYa(6:6:end-6,:);
    subLYb(2:elem,:) = intLYb(6:6:end-6,:);
    subN1a(2:elem,:) = intN1a(6:6:end-6,:);
    subN1b(2:elem,:) = intN1b(6:6:end-6,:);
    subN2a(2:elem,:) = intN2a(6:6:end-6,:);
    subN2b(2:elem,:) = intN2b(6:6:end-6,:);
    
    %compute integral
    PX(:,1:elem) = (intLXa(6:6:end,:)-subLXa)';
    PX(:,2:elem+1) = PX(:,2:elem+1) + intLXb(6:6:end,:)'-subLXb';
    PY(:,1:elem) = intLYa(6:6:end,:)'-subLYa';
    PY(:,2:elem+1) = PY(:,2:elem+1) + intLYb(6:6:end,:)'-subLYb';
    
    PI1(:,1:elem) = (intN1a(6:6:end,:)-subN1a)';
    PI1(:,2:elem+1) = PI1(:,2:elem+1) + (intN1b(6:6:end,:)-subN1b)';
    PI2(:,1:elem) = (intN2a(6:6:end,:)-subN2a)';
    PI2(:,2:elem+1) = PI2(:,2:elem+1) + (intN2b(6:6:end,:)-subN2b)';

end