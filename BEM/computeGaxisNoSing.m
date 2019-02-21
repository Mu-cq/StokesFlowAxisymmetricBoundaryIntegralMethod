%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [G1,G2,Gx,Gy] = computeGaxisNoSing(x,y,x0,y0)

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x0';
    Y0 = y0';
    
    %moltiplicate middle points coordinates
    xMiddle = (x(1:end-1)+x(2:end))/2;  yMiddle = (y(1:end-1)+y(2:end))/2;
    tempx = repmat(xMiddle,6,1);
    XX = reshape(tempx,1,6*elem);
    tempy = repmat(yMiddle,6,1);
    YY = reshape(tempy,1,6*elem);
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX = repmat(GP,1,elem).*reshape((repmat(deltaX/2,6,1)),1,6*elem);
    GPY = repmat(GP,1,elem).*reshape((repmat(deltaY/2,6,1)),1,6*elem);
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,6*elem,1);
    Y0matr = repmat(Y0,6*elem,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %compute gree function and the gradient
    [PHI,PHIX,PHIY] = axisKernelLaplace(globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %multiply times the radius
    PHI = PHI.*globalYmatr;
    PHIX = PHIX.*globalYmatr;
    PHIY = PHIY.*globalYmatr;
    
    temp = reshape(repmat(deltaL/2,6,1),1,6*elem);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',elem,N);
    
    %INTEGRATION
    intP = cumsum(PHI.*manyGW.*manyDelta);
    intPx = cumsum(PHIX.*manyGW.*manyDelta);
    intPy = cumsum(PHIY.*manyGW.*manyDelta);
    
    subP(2:elem,:) = intP(6:6:end-6,:);
    subPx(2:elem,:) = intPx(6:6:end-6,:);
    subPy(2:elem,:) = intPy(6:6:end-6,:);
    
    G1 = intP(6:6:end,:)'-subP';
    Gx = intPx(6:6:end,:)'-subPx';
    Gy = intPy(6:6:end,:)'-subPy';
    
    %compute normal vectors
    no = sqrt(diff(x).^2+diff(y).^2);
    nx = diff(y)./no;
    ny = -diff(x)./no;
    nnx = repmat(nx,N,1);
    nny = repmat(ny,N,1);
    
    %compute kernel for double layer
    G2 = Gx.*nnx + Gy.*nny;

end