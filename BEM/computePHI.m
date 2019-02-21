%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [G,Gx,Gy] = computePHI(x,y)

    %tic

    %number of singularities
    N = numel(x)-1;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = (x(1:end-1)+x(2:end))/2;
    Y0 = (y(1:end-1)+y(2:end))/2;
    
    %moltiplicate sing points coordinates
    tempx = repmat(X0,6,1);
    XX = reshape(tempx,1,6*N);
    tempy = repmat(Y0,6,1);
    YY = reshape(tempy,1,6*N);
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX = repmat(GP,1,N).*reshape((repmat(deltaX/2,6,1)),1,6*N);
    GPY = repmat(GP,1,N).*reshape((repmat(deltaY/2,6,1)),1,6*N);
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0,6*N,1);
    Y0matr = repmat(Y0,6*N,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %CAREFUL BECAUSE I DEFINE GREEN FUNCTION WITH OTHER SIGN
    [P,Px,Py] = gf2D_laplace_fs(globalXmatr,globalYmatr,X0matr,Y0matr,-1);
    
    %singular treatment
    [P,Px,Py,P_an] = poisson2D_ST_const(P,Px,Py,deltaL);

    temp = reshape(repmat(deltaL/2,6,1),1,6*N);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',N,N);
    
    %INTEGRATION
    intP = cumsum(P.*manyGW.*manyDelta);
    intPx = cumsum(Px.*manyGW.*manyDelta);
    intPy = cumsum(Py.*manyGW.*manyDelta);
    
    subP(2:N,:) = intP(6:6:end-6,:);
    subPx(2:N,:) = intPx(6:6:end-6,:);
    subPy(2:N,:) = intPy(6:6:end-6,:);
    
    G = intP(6:6:end,:)'-subP' + P_an;
    Gx = intPx(6:6:end,:)'-subPx';
    Gy = intPy(6:6:end,:)'-subPy';

end