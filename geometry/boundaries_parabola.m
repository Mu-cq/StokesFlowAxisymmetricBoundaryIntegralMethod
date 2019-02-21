%draw closed rectangle

function [a,b] = boundaries_parabola(PARAM)

    L = PARAM.L;
    R = PARAM.R;
    n = PARAM.n;
    m1 = PARAM.m1;
    m2 = PARAM.m2;
    j = PARAM.j;

    %outlet
    [a(1:n+1),b(1:n+1)]=drawline2(L,-0.5*R,L,R,n);
    
    %upper wall
    [a(n+1:n+m1+1),b(n+1:n+m1+1)]=drawline2(L,R,0,R,m1);
    
    %inlet
    %[a(n+m1+1:n+m1+j+1),b(n+m1+1:n+m1+j+1)]=drawline2(0,R,0,-R,j);
    [a(n+m1+1:n+m1+j+1),b(n+m1+1:n+m1+j+1)] = drawline_atan(0,R,0,-R,j);
    
    %lower wall
    [a(n+m1+j+1:n+m1+j+m2+1),b(n+m1+j+1:n+m1+j+m2+1)]=draw_parabola(0,-R,L,-0.5*R,m2);
    
    figure
    plot(a,b,'o-')
    axis equal
    xlabel('x')
    ylabel('y')
    grid on
    
end

