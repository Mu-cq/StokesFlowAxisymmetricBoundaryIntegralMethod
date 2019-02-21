%draw closed rectangle

function [a,b] = boundaries_2Dopen(PARAM)

    L = PARAM.L;
    R = PARAM.R;
    m = PARAM.m;
    
    %upper wall
    [a(1:m+1),b(1:m+1)]=drawline2(L,R,0,R,m);
    
    %lower wall
    [a(m+2:m+m+2),b(m+2:m+m+2)]=drawline2(0,-R,L,-R,m);
    
%     figure
%     plot(a,b,'o-')
%     axis equal
%     xlabel('x')
%     ylabel('y')
    
end

