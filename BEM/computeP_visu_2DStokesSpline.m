%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,PIXX,PIXY,PIYX,PIYY] = computeP_visu_2DStokesSpline(x,y,X0,Y0)

    %tic

    %number of singularities
    N = numel(X0);
    %number of elements
    el = numel(x)-1;
    
    %compute splines coefficients
    [ax, bx, cx, dx, ay, by, cy, dy] = my_spline_periodic (x',y');
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    %many splines coeff
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
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    %h0 = sqrt(bx.*bx+by.*by);
    %h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    manyh = repmat(h,N,1)';
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    X0matr = repmat(X0',6*el,1);
    Y0matr = repmat(Y0',6*el,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);
        
    %not clear what is X0!!!
    [Px,Py,PIxx,PIxy,PIyx,PIyy] = PPI_2DStokes (globalXmatr,globalYmatr,X0matr,Y0matr);
    
    %metrics term
    Px = manyh.*Px;
    Py = manyh.*Py;
    PIxx = manyh.*PIxx;
    PIxy = manyh.*PIxy;
    PIyx = manyh.*PIyx;
    PIyy = manyh.*PIyy;
    
    %INTEGRATION
    %moltiplicate normals
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
    intPIxxa = cumsum(PIxx.*manyGW.*PHIA.*R1);
    intPIxxb = cumsum(PIxx.*manyGW.*PHIB.*R1);
    intPIxya = cumsum(PIxy.*manyGW.*PHIA.*R2);
    intPIxyb = cumsum(PIxy.*manyGW.*PHIB.*R2);
    intPIyxa = cumsum(PIyx.*manyGW.*PHIA.*R1);
    intPIyxb = cumsum(PIyx.*manyGW.*PHIB.*R1);
    intPIyya = cumsum(PIyy.*manyGW.*PHIA.*R2);
    intPIyyb = cumsum(PIyy.*manyGW.*PHIB.*R2);
    
    subPxa(2:el,:) = intPxa(6:6:end-6,:);
    subPxb(2:el,:) = intPxb(6:6:end-6,:);
    subPya(2:el,:) = intPya(6:6:end-6,:);
    subPyb(2:el,:) = intPyb(6:6:end-6,:);
    subPIxxa(2:el,:) = intPIxxa(6:6:end-6,:);
    subPIxxb(2:el,:) = intPIxxb(6:6:end-6,:);
    subPIxya(2:el,:) = intPIxya(6:6:end-6,:);
    subPIxyb(2:el,:) = intPIxyb(6:6:end-6,:);
    subPIyxa(2:el,:) = intPIyxa(6:6:end-6,:);
    subPIyxb(2:el,:) = intPIyxb(6:6:end-6,:);
    subPIyya(2:el,:) = intPIyya(6:6:end-6,:);
    subPIyyb(2:el,:) = intPIyyb(6:6:end-6,:);
    
    %INTEGRATION
    tempB = (intPxb(6:6:end,:)-subPxb)';
    PX = (intPxa(6:6:end,:)-subPxa)' + [tempB(:,end) tempB(:,1:end-1)];
    tempB = (intPyb(6:6:end,:)-subPyb)';
    PY = (intPya(6:6:end,:)-subPya)' + [tempB(:,end) tempB(:,1:end-1)];
    tempB = (intPIxxb(6:6:end,:)-subPIxxb)';
    PIXX = (intPIxxa(6:6:end,:)-subPIxxa)' + [tempB(:,end) tempB(:,1:end-1)];
    tempB = (intPIxyb(6:6:end,:)-subPIxyb)';
    PIXY = (intPIxya(6:6:end,:)-subPIxya)' + [tempB(:,end) tempB(:,1:end-1)];
    tempB = (intPIyxb(6:6:end,:)-subPIyxb)';
    PIYX = (intPIyxa(6:6:end,:)-subPIyxa)' + [tempB(:,end) tempB(:,1:end-1)];
    tempB = (intPIyyb(6:6:end,:)-subPIyyb)';
    PIYY = (intPIyya(6:6:end,:)-subPIyya)' + [tempB(:,end) tempB(:,1:end-1)];
    
    %T = toc

end