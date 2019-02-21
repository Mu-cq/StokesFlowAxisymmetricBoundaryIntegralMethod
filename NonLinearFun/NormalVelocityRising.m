%compute velocity normal to the interface for rising droplet

function uDropFrame = NormalVelocityRising(theta,r,q,lambda,capillary,PARAM)

    %cartesian coordiantes
    a = r'.*cos(theta);
    b = r'.*sin(theta);
    
    %compute solution
    [y,~,~,~,N] = bem(a,b,q+1,lambda,1,capillary,PARAM.way_curv,PARAM.lealpoz,PARAM.elem,PARAM);
    
    %normal velocity IN DROPLET FRAME OF REFERENCE
    ux = y(1:2:end-1);  uy = y(2:2:end);
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    
    %compute velocity of the droplet
    Vdrop = DropVelocityAxis(a,b,u);
    
    %norrmal velocity in droplet frame
    uDropFrame = N(1,:)'.*(ux-Vdrop) + N(2,:)'.*uy;
    
end