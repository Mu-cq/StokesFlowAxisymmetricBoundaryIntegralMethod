%compute normal vector

function [nx,ny] = computeNormalVector(x,y,orderVariable,orderGeometry,SPlinesType)

if orderVariable==0 && orderGeometry==0         % in the midlle of the straight element
   
    no = sqrt(diff(x).^2+diff(y).^2);
    nx = diff(y)./no;
    ny = -diff(x)./no;
    
elseif orderVariable==0 && orderGeometry==1     % in the midlle of the curved element
    
    %compute spline coeff
    if SPlinesType==1
         [~,bx,cx,dx,~,by,cy,dy] = spline_natural(x,y);
    elseif SPlinesType==2
         [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(x,y);
    end
    
    %first and second derivative with splines
    xp = derSplinesConst(bx,cx,dx);
    yp = derSplinesConst(by,cy,dy);
    
    %compute normal vector
    h = sqrt(xp.^2+yp.^2);
    nx = yp./h;
    ny = -xp./h;
    
elseif orderVariable==1 && (orderGeometry==0 || orderGeometry==1)     % in the node (with splines)
    
    %compute spline coeff
    if SPlinesType==1
          [~,bx,cx,dx,~,by,cy,dy] = spline_natural(x,y);
    elseif SPlinesType==2
          [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(x,y);
    end
    
    %first and second derivative with splines
    xp = derSplines(bx,cx,dx);
    yp = derSplines(by,cy,dy);
    
    %compute normal vector
    h = sqrt(xp.^2+yp.^2);
    nx = yp./h;
    ny = -xp./h;
    
end