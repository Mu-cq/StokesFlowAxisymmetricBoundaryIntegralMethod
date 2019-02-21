%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,PI1,PI2] = computeP_visu_const_StokesAX(x,y,X0,Y0)

    %tic

    %number of fields points
    N = numel(X0);
    %number of elements
    el = numel(x)-1;
    
    %normal vector
    noWall = sqrt(diff(x).^2+diff(y).^2);
    nx = -diff(y)./noWall;
    ny = diff(x)./noWall;
    nx = repmat(nx,numel(X0),1);
    ny = repmat(ny,numel(X0),1);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % points where I perform gauss integration
    deltaX = x(2:end)-x(1:end-1);
    deltaY = y(2:end)-y(1:end-1);
    deltaL = sqrt(deltaX.*deltaX+deltaY.*deltaY);
    
    %every Gauss point
    GPX = repmat(GP,1,el).*reshape((repmat(deltaX/2,6,1)),1,6*el);
    GPY = repmat(GP,1,el).*reshape((repmat(deltaY/2,6,1)),1,6*el);
    
    XX = repmat((x(2:end)+x(1:end-1))/2,6,1);
    YY = repmat((y(2:end)+y(1:end-1))/2,6,1);
    globalX = XX(:)'+GPX;
    globalY = YY(:)'+GPY;
    
%     figure
%     plot(globalX,globalY,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
     
    X0matr = repmat(X0,1,6*el)';
    Y0matr = repmat(Y0,1,6*el)';
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %not clear what is X0!!!
    %[Px,Py] = sgf_ax_fs_pressure (globalXmatr,globalYmatr,X0matr,Y0matr);
    [Px,Py,PIxx,PIxy,PIyx,PIyy] = kernelPressureAxis(globalXmatr,globalYmatr,X0matr,Y0matr);

    temp = reshape(repmat(deltaL/2,6,1),6*el,1);
    manyDelta = repmat(temp,1,N);
    manyGW = repmat(GW',el,N);
    
    %INTEGRATION
    intPx = cumsum(Px.*manyGW.*manyDelta);
    intPy = cumsum(Py.*manyGW.*manyDelta);
    intPIxx = cumsum(PIxx.*manyGW.*manyDelta);
    intPIxy = cumsum(PIxy.*manyGW.*manyDelta);
    intPIyx = cumsum(PIyx.*manyGW.*manyDelta);
    intPIyy = cumsum(PIyy.*manyGW.*manyDelta);
    
    subPx(2:el,:) = intPx(6:6:end-6,:);
    subPy(2:el,:) = intPy(6:6:end-6,:);
    subPIxx(2:el,:) = intPIxx(6:6:end-6,:);
    subPIxy(2:el,:) = intPIxy(6:6:end-6,:);
    subPIyx(2:el,:) = intPIyx(6:6:end-6,:);
    subPIyy(2:el,:) = intPIyy(6:6:end-6,:);
    
    PX = intPx(6:6:end,:)'-subPx';
    PY = intPy(6:6:end,:)'-subPy';
    PIXX = intPIxx(6:6:end,:)'-subPIxx';
    PIXY = intPIxy(6:6:end,:)'-subPIxy';
    PIYX = intPIyx(6:6:end,:)'-subPIyx';
    PIYY = intPIyy(6:6:end,:)'-subPIyy';
    
    %sum with normal
    PI1 = PIXX.*nx + PIXY.*ny;
    PI2 = PIYX.*nx + PIYY.*ny;
    
    %T = toc

end