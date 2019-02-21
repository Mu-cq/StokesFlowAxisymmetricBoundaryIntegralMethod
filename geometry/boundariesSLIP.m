%create two vectors conteining the boundaries nodes (first point and last point of the element) once you have given the radius R, the
%lenght L of the pipe, the number of element in each region (n for the
%outlet, m for the wall, j for the inlet).

function [a,b]=boundariesSLIP(R,L,n,m,j,k,delta)

     %h = L*delta;   %length of the slip region
     hh = L*(1-delta);   %length of the no-slip region

     %build the panel at the outlet
     [a(1:n+1),b(1:n+1)]=drawline2(L,0,L,R,n);

     %build the panel at the wall SLIP
     [a(n+1:k+n+1),b(n+1:k+n+1)]=drawline2(L,R,hh,R,k);
     
     %build the panel at the wall NO-SLIP
     [a(n+k+1:m+k+n+1),b(n+k+1:m+k+n+1)]=drawline2(hh,R,0,R,m);

     %build the panel at the inlet
     [a(n+m+k+1:n+m+j+k+1),b(n+m+k+1:n+m+j+k+1)]=drawline2(0,R,0,0,j);
     
%      figure
%      plot(a,b,'o-')
%      axis equal
%      xlabel('x')
%      ylabel('r')
%      title('domain')
%      hold on
end