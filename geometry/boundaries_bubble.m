%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_bubble(PARAM)

    n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
    R_bubble = PARAM.alpha*PARAM.R;

     if R_bubble >= PARAM.R
         b = 0.9*PARAM.R;
         a = R_bubble^3/b^2;
         D = (a-b)/(a+b);
     else
         D = 0;
     end         

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline2(PARAM.L,0,PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),PARAM.n);

     %build the panel at the wall
     [a(n+1:m+n+1),b(n+1:m+n+1)]=drawline2(PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),0,PARAM.R-PARAM.L/2*sin(PARAM.theta),PARAM.m);

     %build the panel at the inlet
     [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline2(0,PARAM.R-PARAM.L/2*sin(PARAM.theta),0,0,PARAM.j);
     
     %build the bubble interface at the first iteration
     [a(n+m+j+2:n+m+j+q+2),b(n+m+j+2:n+m+j+q+2)] = draw_bubble(PARAM.start,0,PARAM.q,D,R_bubble);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end