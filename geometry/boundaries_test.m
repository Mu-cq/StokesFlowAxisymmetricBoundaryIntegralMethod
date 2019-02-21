%draw closed rectangle

function [a,b] = boundaries_test(PARAM)

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
    
    %inner boundary
    [a(n+m+j+m+2:n+2*m+j+2+q),b(n+m+j+m+2:n+2*m+j+2+q)] = draw_test(PARAM.L/2,0,R*PARAM.alpha,q);
    
%     figure
%     plot(a,b,'o-')
%     axis equal
%     xlabel('x')
%     ylabel('y')
    
end

