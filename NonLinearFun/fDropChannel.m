%droplet in channel, normal velocity in drop frame of reference

function vel = fDropChannel(t,xyDrop,xWall,yWall,PARAM)

    display(['t=' num2str(t)])

    %define surface tension
    U = PARAM.Q/pi/PARAM.R^2;
    PARAM.gamma = abs(U)*PARAM.visc2/PARAM.Ca;

    %elements
    n = PARAM.n;    m = PARAM.m;    q = PARAM.q;    j = PARAM.j;
    
    %initialize
    vel = zeros(2*q+2,1);
    
    %general coordinates
    xDrop = xyDrop(1:2:end-1);  yDrop = xyDrop(2:2:end);
    a = [xWall xDrop'];  b = [yWall yDrop'];
    
    %compute velocities
    [y,~,~,~,N] = BEMdropChannelOde(a,b,PARAM);
    
    %drop coordinates and velocities
    aDrop = a(n+m+j+2:end);     bDrop = b(n+m+j+2:end);
    ux = y(2*(n+m+j)+1:2:end-1);    uy = y(2*(n+m+j)+2:2:end);
    
    %normal velocity in drop frame
    uNormal = ux.*N(1,:)' + uy.*N(2,:)';
    Udrop = DropVelocityAxis(aDrop,bDrop,uNormal);
    uNormal = (ux-Udrop).*N(1,:)' + uy.*N(2,:)';
    
    %x and y componenets of normal velocity
    uNormalX = uNormal.*N(1,:)';
    uNormalY = uNormal.*N(2,:)';
    
    %output vector
    vel(1:2:end-1) = uNormalX;
    vel(2:2:end) = uNormalY;

end