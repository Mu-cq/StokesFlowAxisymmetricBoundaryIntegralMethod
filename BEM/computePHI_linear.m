%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [G,G2] = computePHI_linear(x,y,r)

    %tic
    
    %number of singularities
    N = numel(x)-1;
    el = N;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x(1:end-1);
    Y0 = y(1:end-1);
    
    %moltiplicate sing points coordinates
    tempx = repmat(X0,6,1);
    XX = reshape(tempx,1,6*N);
    tempy = repmat(Y0,6,1);
    YY = reshape(tempy,1,6*N);
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %modify GP
    GP = GP+1;
    
    %every Gauss point
    GPX = repmat(GP,1,el).*reshape((repmat(deltaX/2,6,1)),1,6*el);
    GPY = repmat(GP,1,el).*reshape((repmat(deltaY/2,6,1)),1,6*el);
    
    phia = 1-GP/2;
    phib = GP/2;
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %CAREFUL BECAUSE I DEFINE GREEN FUNCTION WITH OTHER SIGN
    [P,Px,Py] = gf2D_laplace_fs(globalXmatr,globalYmatr,X0matr,Y0matr,-1);
    
    %singular treatment
    P_anLin = zeros(numel(X0));
    %[P,Px,Py,P_an] = poisson2D_ST_linear(P,Px,Py,deltaL);


    %INTEGRATION
    temp = reshape(repmat(deltaL/2,6,1),1,6*el);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',el,N);
    PHIA = repmat(phia',el,N);
    PHIB = repmat(phib',el,N);
            
    intPa = cumsum(P.*manyGW.*manyDelta.*PHIA);
    intPb = cumsum(P.*manyGW.*manyDelta.*PHIB);
    intPxa = cumsum(Px.*manyGW.*manyDelta.*PHIA);
    intPxb = cumsum(Px.*manyGW.*manyDelta.*PHIB);
    intPya = cumsum(Py.*manyGW.*manyDelta.*PHIA);
    intPyb = cumsum(Py.*manyGW.*manyDelta.*PHIB);
    
    subPa(2:el,:) = intPa(6:6:end-6,:);
    subPb(2:el,:) = intPb(6:6:end-6,:);
    subPxa(2:el,:) = intPxa(6:6:end-6,:);
    subPxb(2:el,:) = intPxb(6:6:end-6,:);
    subPya(2:el,:) = intPya(6:6:end-6,:);
    subPyb(2:el,:) = intPyb(6:6:end-6,:);
    
    %compute integral (TAKE IN ACCOUNT THAT LAST ELEMENT TOUCHES THE FIRST NODE)    
    temp = (intPb(6:6:end,:)-subPb)';
    G = (intPa(6:6:end,:)-subPa)' + [temp(:,end) temp(:,1:end-1)] + P_anLin;
    
    %multiplicate normals
    R1 = repmat(r(1,:),N,1);
    R2 = repmat(r(2,:),N,1);
    
    temp = (intPxb(6:6:end,:)'-subPxb').*R1 + (intPyb(6:6:end,:)'-subPyb').*R2;
    G2 = (intPxa(6:6:end,:)'-subPxa').*R1 + (intPya(6:6:end,:)'-subPya').*R2+...
        [temp(:,end) temp(:,1:end-1)];
    
    %T = toc

end