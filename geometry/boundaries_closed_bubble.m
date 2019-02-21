%draw closed rectangle

function [a,b] = boundaries_closed_bubble(PARAM)

    L = PARAM.L;
    R = PARAM.R;
    n = PARAM.n;
    m = PARAM.m;
    j = PARAM.j;
    q = PARAM.q;

    %outlet
    [a(1:n+1),b(1:n+1)]=drawline2(L,-R,L,R,n);
    
    %upper wall
    [a(n+1:n+m+1),b(n+1:n+m+1)]=drawline2(L,R,0,R,m);
    
    %inlet
    [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline2(0,R,0,-R,j);
    
    %lower wall
    [a(n+m+j+1:n+m+j+m+1),b(n+m+j+1:n+m+j+m+1)]=drawline2(0,-R,L,-R,m);
    
    %bubble
    [a(n+m+j+m+2:n+2*m+j+2+q),b(n+m+j+m+2:n+2*m+j+2+q)] = draw_closed_circle(PARAM.x0,0,PARAM.alpha*R,q,0,0);
    
    %modify bubble shape if the channle is bigger than the channel
    if PARAM.alpha>1
    h = R*(0.643*(3*PARAM.Ca)^(2/3)/(1+2.5*0.643*(3*PARAM.Ca)^(2/3)));   %from Bretherton
    if PARAM.visc >= 0.1
        h = h*4^(2/3);  %Odge correction for viscosity
    end
    A = pi*PARAM.alpha^2*R^2;
    l = (A-pi*(R-h)^2)/2/(R-h); %width inner part
    
    %choose bumber of elements per part
    cut = l/(pi*(R-h)+l);
    q_str = ceil(cut*q/2);    %elements for straigth part
    q_round = q/2-q_str; %elements for rounded part
    
    %divide each part
    q_right = ceil(q_round/2);  q_left = q_round-q_right;
    q_up = ceil(q_str);
    
    %build bubble
    theta = 0:pi/q_right/2:pi/2;
    Xright = (R-h)*cos(theta)+L/2+l/2;   Yright = (R-h)*sin(theta);
    theta = pi/2:pi/q_left/2:pi;
    Xleft = (R-h)*cos(theta)+L/2-l/2;   Yleft = (R-h)*sin(theta);
    
    Xup = (l/2:-l/q_up:-l/2)+L/2; Yup = (R-h)*ones(1,numel(Xup));
    
    %create upper part
    aUp = [Xright(1:end-1) Xup(1:end-1) Xleft];
    bUp = [Yright(1:end-1) Yup(1:end-1) Yleft];
    
    %lower part is created by mirroring
    a(n+m+j+m+2:n+2*m+j+2+q) = [aUp(1:end-1) flip(aUp)];
    b(n+m+j+m+2:n+2*m+j+2+q) = [bUp(1:end-1) -flip(bUp)];
    
    
    end
    
%     figure
%     plot(a,b,'o-')
%     axis equal
%     xlabel('x')
%     ylabel('y')
    
end

