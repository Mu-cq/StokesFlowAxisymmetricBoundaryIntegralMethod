%create two vectors conteining the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundaries_smooth(R,L,n,m,j)

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline2(L,0,L,R,n);

     %build the panel at the wall
     [a(n+1:m+n+1),b(n+1:m+n+1)]=drawline2(L,R,0,R,m);

     %build the panel at the inlet
     [a(n+m+1:n+m+j+1),b(n+m+1:n+m+j+1)]=drawline2(0,R,0,0,j);
     
     %smooth corners
     a(n+1) = (a(n+2)+a(n))/2;
     b(n+1) = (b(n+2)+b(n))/2;
     
     %smooth corners
     a(n+m+1) = (a(n+m+2)+a(n+m))/2;
     b(n+m+1) = (b(n+m+2)+b(n+m))/2;
     
%      figure
%      plot(a,b,'o-')
%      axis equal
%      xlabel('x')
%      ylabel('r')
%      title('domain')
%      hold on
end