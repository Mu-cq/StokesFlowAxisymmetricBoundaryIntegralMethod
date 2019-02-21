%compute velovity normal to the interface

function Un = DropNormalVelocity(x,y,u,v)

    %compute splines coeff
    [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(x,y);
    
    %compute normal vector
    [nx,ny] = DropNormalVector(bx,cx,dx,by,cy,dy);
    
    %compute normal velocity
    Un = nx.*u + ny.*v;

end