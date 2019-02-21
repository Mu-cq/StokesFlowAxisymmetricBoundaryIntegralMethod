%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,PI1,PI2] = computeP_AX_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0,Y0)

    %tic

    %number of singularities
    N = numel(X0);
    %number of element
    el = numel(ax);
    
    %allocation
    PX = zeros(N,el+1);
    PY = zeros(N,el+1);
    PI1 = zeros(N,el+1);
    PI2 = zeros(N,el+1);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %spline confficient for vectorial
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
    
    %hat functions
    phia = 1-GPt;
    phib = GPt;
    
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
    %h0 = sqrt(bx.*bx+by.*by);
    %h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
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

    [Px,Py,Pixx,Pixy,~,Piyy] = kernelPressureAxis(globalXmatr,globalYmatr,X0matr,Y0matr);
            
    Px = manyh.*Px;
    Py = manyh.*Py;
    Pixx = manyh.*Pixx;
    Pixy = manyh.*Pixy;
    Piyx = Pixy;
    Piyy = manyh.*Piyy;
       
    %moltiplicate normals for double layer
    R1 = repmat(nx,N,1)';
    R2 = repmat(ny,N,1)';
    
    %for integration
    manyGW = repmat(GW',el,N);
    PHIA = repmat(phia',el,N);
    PHIB = repmat(phib',el,N);
    
    intPxa = cumsum(Px.*manyGW.*PHIA);
    intPxb = cumsum(Px.*manyGW.*PHIB);
    intPya = cumsum(Py.*manyGW.*PHIA);
    intPyb = cumsum(Py.*manyGW.*PHIB);
    intPixxa = cumsum(Pixx.*manyGW.*PHIA.*R1);
    intPixxb = cumsum(Pixx.*manyGW.*PHIB.*R1);
    intPixya = cumsum(Pixy.*manyGW.*PHIA.*R2);
    intPixyb = cumsum(Pixy.*manyGW.*PHIB.*R2);
    intPiyxa = cumsum(Piyx.*manyGW.*PHIA.*R1);
    intPiyxb = cumsum(Piyx.*manyGW.*PHIB.*R1);
    intPiyya = cumsum(Piyy.*manyGW.*PHIA.*R2);
    intPiyyb = cumsum(Piyy.*manyGW.*PHIB.*R2);
    
    subPxa(2:el,:) = intPxa(6:6:end-6,:);
    subPxb(2:el,:) = intPxb(6:6:end-6,:);
    subPya(2:el,:) = intPya(6:6:end-6,:);
    subPyb(2:el,:) = intPyb(6:6:end-6,:);
    subPixxa(2:el,:) = intPixxa(6:6:end-6,:);
    subPixxb(2:el,:) = intPixxb(6:6:end-6,:);
    subPixya(2:el,:) = intPixya(6:6:end-6,:);
    subPixyb(2:el,:) = intPixyb(6:6:end-6,:);
    subPiyxa(2:el,:) = intPiyxa(6:6:end-6,:);
    subPiyxb(2:el,:) = intPiyxb(6:6:end-6,:);
    subPiyya(2:el,:) = intPiyya(6:6:end-6,:);
    subPiyyb(2:el,:) = intPiyyb(6:6:end-6,:);
    
    %compute integral single layer integral
    PX(:,1:el) = (intPxa(6:6:end,:)-subPxa)';
    PX(:,2:el+1) = PX(:,2:el+1) + (intPxb(6:6:end,:)-subPxb)';
    PY(:,1:el) = (intPya(6:6:end,:)-subPya)';
    PY(:,2:el+1) = PY(:,2:el+1) + (intPyb(6:6:end,:)-subPyb)';
    
    %compute double layer integral
    PI1(:,1:el) = (intPixxa(6:6:end,:)-subPixxa + intPixya(6:6:end,:)-subPixya)';
    PI1(:,2:el+1) = PI1(:,2:el+1) + (intPixxb(6:6:end,:)-subPixxb + intPixyb(6:6:end,:)-subPixyb)';
    PI2(:,1:el) = (intPiyxa(6:6:end,:)-subPiyxa + intPiyya(6:6:end,:)-subPiyya)';
    PI2(:,2:el+1) = PI1(:,2:el+1) + (intPiyxb(6:6:end,:)-subPiyxb + intPiyyb(6:6:end,:)-subPiyyb)'; 
    
    %T = toc

end