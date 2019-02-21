%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_narrow_large2(PARAM)

     n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;   p = PARAM.p;

     if PARAM.alpha*PARAM.R >= PARAM.R
         b = 0.9*R;
         a = R_bubble^3/b^2;
         D = (a-b)/(a+b);
     else
         D = 0;
     end

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline2(PARAM.L-PARAM.R*(PARAM.large-1)*2,0,PARAM.L-PARAM.R*(PARAM.large-1)*2,PARAM.R,n);

     %build the panel at the wall
     [a1,b1]=drawline2(PARAM.L-PARAM.R*(PARAM.large-1)*2,PARAM.R,...
         PARAM.L-PARAM.R*(PARAM.large-1)*2-PARAM.leng,PARAM.R,m/PARAM.L*PARAM.leng);
     [a2,b2]=drawline2(PARAM.L-PARAM.R*(PARAM.large-1)*2-PARAM.leng,PARAM.R,...
         PARAM.L-PARAM.R*(PARAM.large-1)*2-PARAM.leng,PARAM.R*PARAM.large,m*PARAM.R*(PARAM.large-1)/PARAM.L);
     [a3,b3]=drawline2(PARAM.L-PARAM.R*(PARAM.large-1)*2-PARAM.leng,PARAM.R*PARAM.large,...
         PARAM.L-PARAM.R*(PARAM.large-1)*2-2*PARAM.leng,PARAM.R*PARAM.large,m/PARAM.L*PARAM.leng);
     [a4,b4]=drawline2(PARAM.L-PARAM.R*(PARAM.large-1)*2-2*PARAM.leng,PARAM.R*PARAM.large,...
         PARAM.L-PARAM.R*(PARAM.large-1)*2-2*PARAM.leng,PARAM.R,m*PARAM.R*(PARAM.large-1)/PARAM.L);
     [a5,b5]=drawline2(PARAM.L-PARAM.R*(PARAM.large-1)*2-2*PARAM.leng,PARAM.R,...
         PARAM.L-PARAM.R*(PARAM.large-1)*2-3*PARAM.leng,PARAM.R,m/PARAM.L*PARAM.leng);
     a(n+1:m+n+1) = [a1 a2(2:end) a3(2:end) a4(2:end) a5(2:end)];
     b(n+1:m+n+1) = [b1 b2(2:end) b3(2:end) b4(2:end) b5(2:end)];

     %build the panel at the inlet
     [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline2(0,PARAM.R,0,0,j);
     
     %build first bubble interface at the first iteration
     [a(n+m+j+2:n+m+j+q+2),b(n+m+j+2:n+m+j+q+2)] = draw_bubble(1,0,q,D,PARAM.R*PARAM.alpha);
     
     %build second bubble interface at the first iteration
     [a(n+m+j+q+3:n+m+j+q+p+3),b(n+m+j+q+3:n+m+j+q+p+3)] = draw_bubble(1+PARAM.dist,0,p,D,PARAM.R*PARAM.alpha);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end