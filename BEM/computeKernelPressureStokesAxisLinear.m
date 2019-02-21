%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,DL1,DL2] = computeKernelPressureStokesAxisLinear(x,y,x0,y0)

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
    %initialize
    PX = zeros(N,elem+1);
    PY = zeros(N,elem+1);
    DL1 = zeros(N,elem+1);
    DL2 = zeros(N,elem+1);
    
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
    [LX,LY,NXX,NXY,NYX,NYY] = kernelPressureAxis(globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %mutuply kernnel of doublelayer times the normal
    PI1 = NXX.*nnx + NXY.*nny;
    PI2 = NYX.*nnx + NYY.*nny;
    
    %INTEGRATION
    temp = reshape(repmat(deltaL/2,6,1),1,6*elem);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',elem,N);
    PHIA = repmat(phia',elem,N);
    PHIB = repmat(phib',elem,N);
    
    %INTEGRATION
    intLXa = cumsum(LX.*manyGW.*manyDelta.*PHIA);
    intLXb = cumsum(LX.*manyGW.*manyDelta.*PHIB);
    intLYa = cumsum(LY.*manyGW.*manyDelta.*PHIA);
    intLYb = cumsum(LY.*manyGW.*manyDelta.*PHIB);
    intPI1a = cumsum(PI1.*manyGW.*manyDelta.*PHIA);
    intPI1b = cumsum(PI1.*manyGW.*manyDelta.*PHIB);
    intPI2a = cumsum(PI2.*manyGW.*manyDelta.*PHIA);
    intPI2b = cumsum(PI2.*manyGW.*manyDelta.*PHIB);
    
    subLXa(2:elem,:) = intLXa(6:6:end-6,:);
    subLXb(2:elem,:) = intLXb(6:6:end-6,:);
    subLYa(2:elem,:) = intLYa(6:6:end-6,:);
    subLYb(2:elem,:) = intLYb(6:6:end-6,:);
    subPI1a(2:elem,:) = intPI1a(6:6:end-6,:);
    subPI1b(2:elem,:) = intPI1b(6:6:end-6,:);
    subPI2a(2:elem,:) = intPI2a(6:6:end-6,:);
    subPI2b(2:elem,:) = intPI2b(6:6:end-6,:);
    
    %compute integral
    PX(:,1:elem) = (intLXa(6:6:end,:)-subLXa)';
    PX(:,2:elem+1) = PX(:,2:elem+1) + intLXb(6:6:end,:)'-subLXb';
    PY(:,1:elem) = intLYa(6:6:end,:)'-subLYa';
    PY(:,2:elem+1) = PY(:,2:elem+1) + intLYb(6:6:end,:)'-subLYb';
    
    DL1(:,1:elem) = (intPI1a(6:6:end,:)-subPI1a)';
    DL1(:,2:elem+1) = DL1(:,2:elem+1) + (intPI1b(6:6:end,:)-subPI1b)';
    DL2(:,1:elem) = (intPI2a(6:6:end,:)-subPI2a)';
    DL2(:,2:elem+1) = DL2(:,2:elem+1) + (intPI2b(6:6:end,:)-subPI2b)';
    
end