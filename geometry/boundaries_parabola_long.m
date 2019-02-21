%draw closed rectangle

function [a,b] = boundaries_parabola_long(PARAM)

    L = PARAM.L;
    R = PARAM.R;
    n = PARAM.n;
    m1 = PARAM.m1;
    m2 = PARAM.m2;
    m3 = PARAM.m3;
    j = PARAM.j;

    %outlet
    [a(1:n+1),b(1:n+1)]=drawline2(PARAM.R+L,-0.5*R,R+L,R,n);
    
    %upper wall
    [a(n+1:n+m1+1),b(n+1:n+m1+1)]=drawline2(R+L,R,0,R,m1);
    
    %inlet
    [a(n+m1+1:n+m1+j+1),b(n+m1+1:n+m1+j+1)]=drawline2(0,R,0,-R,j);
    %[a(n+m1+1:n+m1+j+1),b(n+m1+1:n+m1+j+1)] = drawline_atan(0,R,0,-R,j);
    
    %lower wall straigth
    [a(n+m1+j+1:n+m1+m2+j+1),b(n+m1+j+1:n+m1+m2+j+1)]=drawline2(0,-R,R,-R,m2);
    
    %lower wall parabola
    [a(n+m1+m2+j+1:n+m1+m2+m3+j+1),b(n+m1+m2+j+1:n+m1+m2+m3+j+1)]=draw_parabola(R,-R,L+R,-0.5*R,m3);
    
    figure
    plot(a,b,'o-')
    axis equal
    xlabel('x')
    ylabel('y')
    grid on
    
end

