%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [G1,G2] = computeKernelLaplaceAxis2(x,y,x0,y0,indIn,indOut,ST)

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
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
    
    %moltiplicate middle points coordinates
    xMiddle = (x(1:end-1)+x(2:end))/2;  yMiddle = (y(1:end-1)+y(2:end))/2;
    tempx = repmat(xMiddle,6,1);
    XX = reshape(tempx,1,6*elem);
    tempy = repmat(yMiddle,6,1);
    YY = reshape(tempy,1,6*elem);
    
    %lenght of the elements
    deltaX = diff(x);
    deltaY = diff(y);
    deltaL = sqrt(deltaX.^2+deltaY.^2);
    
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
    PHIXnx_PHIYny = PHIX.*globalYmatr.*nnx + PHIY.*globalYmatr.*nny;
    
    %singular treatment
    if ST==1
        range = indIn:indOut;
        [PHI(:,range),PHIXnx_PHIYny(:,range),P_an,PxPy_an] = laplaceST2(PHI(:,range),PHIXnx_PHIYny(:,range),deltaL,globalXmatr(:,range),globalYmatr(:,range),x0(range),y0(range),ny);
        %P_an = P_an.*
    end
    
    temp = reshape(repmat(deltaL/2,6,1),1,6*elem);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',elem,N);
    
    %INTEGRATION
    intP = cumsum(PHI.*manyGW.*manyDelta);
    intPHIXnx_PHIYny = cumsum(PHIXnx_PHIYny.*manyGW.*manyDelta);
    
    subP(2:elem,:) = intP(6:6:end-6,:);
    subPHIXnx_PHIYny(2:elem,:) = intPHIXnx_PHIYny(6:6:end-6,:);
    
    G1 = intP(6:6:end,:)'-subP';
    G2 = intPHIXnx_PHIYny(6:6:end,:)'-subPHIXnx_PHIYny';
    
    %add anaytical integration of the singularity
    if ST==1
        G1(range,:) = G1(range,:) + P_an;
        G2(range,:) = G2(range,:) + PxPy_an;
    end

end