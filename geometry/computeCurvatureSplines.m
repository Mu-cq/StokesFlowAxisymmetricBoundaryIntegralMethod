%compute normal vector

function [k1,k2] = computeCurvatureSplines(x,y,orderVariable)

if orderVariable==0         % in the midlle of the straight element
   
    error('It is not possible to compute the curvature in the middle of the element')
    
elseif orderVariable==1     % in the node (with splines)
    
    %compute spline coeff
    [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(x, y);
    
    %first and second derivative with splines
    [xp,xpp] = derSplines(bx,cx,dx);
    [yp,ypp] = derSplines(by,cy,dy);
    
    %compute normal vector
    h = sqrt(xp.^2+yp.^2);
    ny = -xp./h;
    
    %curvature in the mid plane
    k1 = -(xpp.*yp-ypp.*xp)./h.^3;
    k2 = ny./y;
    k2([1 end]) = k1([1 end]);
    
end