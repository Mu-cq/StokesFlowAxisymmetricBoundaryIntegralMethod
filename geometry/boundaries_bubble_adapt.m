%create two vectors containing the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_bubble_adapt(R,L,n,m,j,q,theta,R_bubble)

     if R_bubble>=R
         b = 0.9*R;
         a = R_bubble^3/b^2;
         D = (a-b)/(a+b);
     else
         D = 0;
     end         

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline2(L,0,L,R+L/2*sin(theta),n);

     %build the panel at the wall
     [a(n+1:m+n+1),b(n+1:m+n+1)]=drawline_adapt(L,R+L/2*sin(theta),0,R-L/2*sin(theta),m);

     %build the panel at the inlet
     [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline2(0,R-L/2*sin(theta),0,0,j);
     
     %build the bubble interface at the first iteration
     [a(n+m+j+2:n+m+j+q+2),b(n+m+j+2:n+m+j+q+2)] = draw_bubble(L/2,0,q,D,R_bubble);
     
%      figure
%      plot(a,b,'o-')
%      xlabel('x')
%      ylabel('r')
%      axis equal
%      title('domain')
%      hold on
end