%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [G1,G2] = computeKernelLaplaceAxisLinear(x,y,x0,y0,indIn,indOut,ST)

    error('Bug here, problem when the singularity appraches the axis: use piecewise consatnt elements')

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
    %initialize
    G1 = zeros(N,elem+1);
    G2 = zeros(N,elem+1);
    
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
    [PHI,PHIX,PHIY] = axisKernelLaplace(globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %multiply times the radius
    PHI = PHI.*globalYmatr;
    PHIXnx_PHIYny = PHIX.*globalYmatr.*nnx + PHIY.*globalYmatr.*nny;
      
    %singular treatment
    if ST==1
        error('not implmented')
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
    intPHIa = cumsum(PHI.*manyGW.*manyDelta.*PHIA);
    intPHIb = cumsum(PHI.*manyGW.*manyDelta.*PHIB);
    intPHIXnx_PHIYnya = cumsum(PHIXnx_PHIYny.*manyGW.*manyDelta.*PHIA);
    intPHIXnx_PHIYnyb = cumsum(PHIXnx_PHIYny.*manyGW.*manyDelta.*PHIB);
    
    subPHIa(2:elem,:) = intPHIa(6:6:end-6,:);
    subPHIb(2:elem,:) = intPHIb(6:6:end-6,:);
    subPHIXnx_PHIYnya(2:elem,:) = intPHIXnx_PHIYnya(6:6:end-6,:);
    subPHIXnx_PHIYnyb(2:elem,:) = intPHIXnx_PHIYnyb(6:6:end-6,:);
    
    %compute integral
    G1(:,1:elem) = (intPHIa(6:6:end,:)-subPHIa)';
    G1(:,2:elem+1) = G1(:,2:elem+1) + intPHIb(6:6:end,:)'-subPHIb';
    G2(:,1:elem) = intPHIXnx_PHIYnya(6:6:end,:)'-subPHIXnx_PHIYnya';
    G2(:,2:elem+1) = G2(:,2:elem+1) + intPHIXnx_PHIYnyb(6:6:end,:)'-subPHIXnx_PHIYnyb';
    
    %add anaytical integration of the singularity
    if ST==1
        
        N = numel(deltaL)+1;
        
        range1 = (2:N-2)+indIn(1)-1;
        range2 = (2:N-2)+indIn(1);
        
        %add analtyical integration
        GXX(range1,range1) = GXX(range1,range1) + anA;
        GXX(range1,range2) = GXX(range1,range2) + anB;
        GXX(range2,range1) = GXX(range2,range1) + anB;
        GXX(range2,range2) = GXX(range2,range2) + anA;
        
        GYY(range1,range1) = GYY(range1,range1) + anA;
        GYY(range1,range2) = GYY(range1,range2) + anB;
        GYY(range2,range1) = GYY(range2,range1) + anB;
        GYY(range2,range2) = GYY(range2,range2) + anA;
        
        %because singularity on the axis don't has to be treated
        GXX(range(2),range(1)) = GXX(range(2),range(1)) - deltaL(1)*log(deltaL(1))+deltaL(1);
        GXX(range(2),range(2)) = GXX(range(2),range(2)) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
        GXX(range(end-1),range(end-1)) = GXX(range(end-1),range(end-1)) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
        GXX(range(end-1),range(end)) = GXX(range(end-1),range(end)) - deltaL(N-1)*log(deltaL(range(end-1)))+deltaL(N-1);
        
        GYY(range(2),range(1)) = GYY(range(2),range(1)) - deltaL(1)*log(deltaL(1))+deltaL(1);
        GYY(range(2),range(2)) = GYY(range(2),range(2)) - deltaL(1)*log(deltaL(1))+2*deltaL(1);
        GYY(range(end-1),range(end-1)) = GYY(range(end-1),range(end-1)) - deltaL(N-1)*log(deltaL(N-1))+2*deltaL(N-1);
        GYY(range(end-1),range(end)) = GYY(range(end-1),range(end)) - deltaL(N-1)*log(deltaL(N-1))+deltaL(N-1);
    
    end

end