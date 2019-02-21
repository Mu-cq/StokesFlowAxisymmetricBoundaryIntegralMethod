%droplet in channel, normal velocity in drop frame of reference

function vel = fDropChannelGhostPoint(t,xyDrop,xWall,yWall,PARAM)

    display(['t=' num2str(t)])

    %define surface tension
    U = PARAM.Q/pi/PARAM.R^2;
    PARAM.gamma = abs(U)*PARAM.visc2/PARAM.Ca;
    
    %elements
    n = PARAM.n;    m = PARAM.m;    q = PARAM.q;    j = PARAM.j;
    
    %initialize
    vel = zeros(2*q,1);
    
    %initial volume
    V0 = 4/3*pi*PARAM.alpha^3;
    
    %general coordinates
    xDrop = xyDrop(1:2:end-1);  yDrop = xyDrop(2:2:end);
    %compute location of ghost point
    fVolume = @(rho) ModifyVolume4(xDrop',yDrop',rho,V0);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    xEnd = fsolve(fVolume,xDrop(end),options);
    
    a = [xWall xDrop' xEnd];  b = [yWall yDrop' 0];
    
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
    vel(1:2:end-1) = uNormalX(1:end-1);
    vel(2:2:end) = uNormalY(1:end-1);
    
    %plot drop
    figure(2)
    plot([xDrop' xEnd],[yDrop' 0],'x-')
    grid on
    axis equal
    axis([min(xDrop)-0.5 max(xDrop)+0.5 0 1])
    drawnow

end