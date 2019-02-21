%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22,T1,T2,D1,D2] = computeGT_isoparametric(x,y,ax,ay,bx,by,cx,cy,dx,dy,PARAM)

    %tic

    %number of nodes
    N = numel(x);
    %number of element
    el = N-1;
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    D1 = zeros(N);
    D2 = zeros(N);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    aax = reshape((repmat(ax,6,1)),1,6*el);
    bbx = reshape((repmat(bx,6,1)),1,6*el);
    ccx = reshape((repmat(cx,6,1)),1,6*el);
    ddx = reshape((repmat(dx,6,1)),1,6*el);
    aay = reshape((repmat(ay,6,1)),1,6*el);
    bby = reshape((repmat(by,6,1)),1,6*el);
    ccy = reshape((repmat(cy,6,1)),1,6*el);
    ddy = reshape((repmat(dy,6,1)),1,6*el);    
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;
    GW = GW/2;
    
    phia = 1-GPt;
    phib = GPt;

    % point where the variable will be stored
    X0 = x;
    Y0 = y;
    
    %points where I perform gauss integration
    GPglobal = repmat(GPt,1,el);
    
    %global coordianted retrived on the splines
    eta = aax+bbx.*GPglobal+ccx.*GPglobal.*GPglobal+ddx.*GPglobal.*GPglobal.*GPglobal;
    beta = aay+bby.*GPglobal+ccy.*GPglobal.*GPglobal+ddy.*GPglobal.*GPglobal.*GPglobal;
    deta = bbx+2*ccx.*GPglobal+3*ddx.*GPglobal.*GPglobal;
    dbeta = bby+2*ccy.*GPglobal+3*ddy.*GPglobal.*GPglobal;
    
    %chabge sign for my convention
    nx = dbeta./sqrt(deta.*deta+dbeta.*dbeta);
    ny = -deta./sqrt(deta.*deta+dbeta.*dbeta);
    
%     figure
%     plot(eta,nx)
%     hold on
%     plot(eta,ny,'r')
%     hold off
    
%     figure
%     plot(eta,beta,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off

%     h = zeros(1,el*6);
%     h0 = zeros(el,1);
%     h1 = zeros(el,1);
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    h0 = sqrt(bx.*bx+by.*by);
    h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    manyh = repmat(h,N,1)';
    
%     figure
%     plot(x(1:end-1),h0,'o-')
%     hold on
%     plot(x(2:end),h1,'r-o')
%     plot(eta,h,'k-o')
%     hold off

%     figure
%     plot(h,'o-')

    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);

    %compute green's functions
    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY,PXX,PXY,PYX,PYY] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    %compute weights for isoparametric
    weightA = PARAM.weightA;
    weightB = PARAM.weightB;
    weightC = PARAM.weightC;
    weightD = PARAM.weightD;
    
    %repeat weights
    weightA = reshape(repmat(weightA(:)',6,1),6*(numel(x)-1),numel(x));
    weightB = reshape(repmat(weightB(:)',6,1),6*(numel(x)-1),numel(x));
    weightC = reshape(repmat(weightC(:)',6,1),6*(numel(x)-1),numel(x));
    weightD = reshape(repmat(weightD(:)',6,1),6*(numel(x)-1),numel(x));
    
    %build matric of weights
    manyGPglobal = repmat(GPglobal',1,numel(x));
    Z = weightA + weightB.*manyGPglobal + weightC.*manyGPglobal.^2 + weightD.*manyGPglobal.^3;
    
%     figure
%     plot(Z)
%     grid on
    
    SXX = manyh.*SXX;
    SXY = manyh.*SXY;
    SYX = manyh.*SYX;
    SYY = manyh.*SYY;
%     QXXX = manyh.*QXXX;
%     QXXY = manyh.*QXXY;
%     QXYX = manyh.*QXYX;
%     QXYY = manyh.*QXYY;
%     QYXX = manyh.*QYXX;
%     QYXY = manyh.*QYXY;
%     QYYX = manyh.*QYYX;
%     QYYY = manyh.*QYYY;
%     PXX = manyh.*PXX;
%     PXY = manyh.*PXY;
%     PYX = manyh.*PYX;
%     PYY = manyh.*PYY;
    
%     figure
%     surf(SYY)
    
    %singularity treatment for linear element REPLACE WITH VECTORIAL
    %PRODUCT %%%%%
%      for i = 2:N-2
% 
%        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i))))';
%        SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i))))';
%               
%        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i)).^2+(beta(1+6*(i-1):6*i)-Y0(i)).^2)./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i))))';
%        SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + (2*log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)).*h(1+6*(i-1):6*i)...
%            - 2*(log(sqrt((eta(1+6*(i-1):6*i)-X0(i+1)).^2+(beta(1+6*(i-1):6*i)-Y0(i+1)).^2)./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i))))';
% 
%     end
%     
%        SXX(1:6,2) = SXX(1:6,2) + (2*log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)).*h(1:6)...
%            - 2*(log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
%        SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) + (2*log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)).*h(1+6*(N-2):6*(N-1))...
%            - 2*(log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1))))';
%     
%        SYY(1:6,2) = SYY(1:6,2) + (2*log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)).*h(1:6)...
%            - 2*(log(sqrt((eta(1:6)-X0(2)).^2+(beta(1:6)-Y0(2)).^2)./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
%        SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) + (2*log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)).*h(1+6*(N-2):6*(N-1))...
%            - 2*(log(sqrt((eta(1+6*(N-2):6*(N-1))-X0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-Y0(N-1)).^2)./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1))))';
       
    %moltiplicate normals
    R1 = repmat(nx,N,1)';
    R2 = repmat(ny,N,1)';
    

    manyGW = repmat(GW',el,1);

    for i = 1:N
        
        for j = 1:N
            
            GXX(j,i) = sum(SXX(:,j).*manyGW.*Z(:,i));
            GXY(j,i) = sum(SXY(:,j).*manyGW.*Z(:,i));
            GYX(j,i) = sum(SYX(:,j).*manyGW.*Z(:,i));
            GYY(j,i) = sum(SYY(:,j).*manyGW.*Z(:,i));
            
%             if i==j
%                 
%             display(i)
%             
%             %check
%             figure(10)
%             plot(globalXmatr(:,1),SXX(:,j));
%             hold on
%             plot(globalXmatr(:,1),Z(:,i));
%             plot(globalXmatr(:,1),SXX(:,j).*Z(:,i));
%             xlabel('x')
%             grid on
%             legend('SXX','Z','product')
%             hold off
%             
%             end
        
        end
        
    end
    
%     A11(:,1:N-1) = (intQXXXa(6:6:end,:)-subQXXXa + intQXXYa(6:6:end,:)-subQXXYa)';
%     A11(:,2:N) = A11(:,2:N) + (intQXXXb(6:6:end,:)-subQXXXb + intQXXYb(6:6:end,:)-subQXXYb)';
%     A12(:,1:N-1) = (intQXYXa(6:6:end,:)-subQXYXa + intQXYYa(6:6:end,:)-subQXYYa)';
%     A12(:,2:N) = A12(:,2:N) + (intQXYXb(6:6:end,:)-subQXYXb + intQXYYb(6:6:end,:)-subQXYYb)';
%     A21(:,1:N-1) = (intQYXXa(6:6:end,:)-subQYXXa + intQYXYa(6:6:end,:)-subQYXYa)';
%     A21(:,2:N) = A21(:,2:N) + (intQYXXb(6:6:end,:)-subQYXXb + intQYXYb(6:6:end,:)-subQYXYb)';
%     A22(:,1:N-1) = (intQYYXa(6:6:end,:)-subQYYXa + intQYYYa(6:6:end,:)-subQYYYa)';
%     A22(:,2:N) = A22(:,2:N) + (intQYYXb(6:6:end,:)-subQYYXb + intQYYYb(6:6:end,:)-subQYYYb)';
%     
%     T1 = A11;
%     T2 = A21;
%     
%     D1(:,1:N-1) = (intPXXa(6:6:end,:)-subPXXa + intPXYa(6:6:end,:)-subPXYa)';
%     D1(:,2:N) = D1(:,2:N) + (intPXXb(6:6:end,:)-subPXXb + intPXYb(6:6:end,:)-subPXYb)';
%     D2(:,1:N-1) = (intPYXa(6:6:end,:)-subPYXa + intPYYa(6:6:end,:)-subPYYa)';
%     D2(:,2:N) = D2(:,2:N) + (intPYXb(6:6:end,:)-subPYXb + intPYYb(6:6:end,:)-subPYYb)';
    
    %singularity treatment
%     GXX(2:N-2,2:N-2) = GXX(2:N-2,2:N-2) + diag(3/2*h0(2:N-2));
%     GXX(2:N-2,3:N-1) = GXX(2:N-2,3:N-1) + diag(1/2*h0(2:N-2));
%     GXX(3:N-1,2:N-2) = GXX(3:N-1,2:N-2) + diag(1/2*h1(2:N-2));
%     GXX(3:N-1,3:N-1) = GXX(3:N-1,3:N-1) + diag(3/2*h1(2:N-2));
%     
%     GYY(2:N-2,2:N-2) = GYY(2:N-2,2:N-2) + diag(3/2*h0(2:N-2));
%     GYY(2:N-2,3:N-1) = GYY(2:N-2,3:N-1) + diag(1/2*h0(2:N-2));
%     GYY(3:N-1,2:N-2) = GYY(3:N-1,2:N-2) + diag(1/2*h1(2:N-2));
%     GYY(3:N-1,3:N-1) = GYY(3:N-1,3:N-1) + diag(3/2*h1(2:N-2));
%     
%     %because singularity on the axis don't has to be treated
%     GXX(2,1) = GXX(2,1) + 1/2*h1(1);
%     GXX(2,2) = GXX(2,2) + 3/2*h1(1);
%     GXX(N-1,N-1) = GXX(N-1,N-1) + 3/2*h0(el);
%     GXX(N-1,N) = GXX(N-1,N) + 1/2*h0(el);
%     
%     GYY(2,1) = GYY(2,1) + 1/2*h1(1);
%     GYY(2,2) = GYY(2,2) + 3/2*h1(1);
%     GYY(N-1,N-1) = GYY(N-1,N-1) + 3/2*h0(el);
%     GYY(N-1,N) = GYY(N-1,N) + 1/2*h0(el);
    
    %T = toc

end