%compute velocity normal to the interface in an extensional flow

function uNorm = NormalVelocityNormRisingXY(XY,q,lambda,capillary,PARAM)

    %cartesian coordiantes
    a = XY([1 2:2:end]);
    b = [0 XY(3:2:end-1) 0];
    
    %compute solution
    [y,~,~,~,N] = bem(a,b,q+1,lambda,1,capillary,PARAM.way_curv,PARAM.lealpoz,PARAM.elem,PARAM);
    
    %normal velocity IN DROPLET FRAME OF REFERENCE
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %compute velocity of the droplet
    Vdrop = DropVelocityAxis(a,b,u);
    
    %norrmal velocity in droplet frame
    uDropFrame = N(1,:)'.*(ux-Vdrop) + N(2,:)'.*uy;
    
    %integral of the velocity on the surface
    uNorm = int_axis_spline_symmetric(a,b,uDropFrame.^2);

end