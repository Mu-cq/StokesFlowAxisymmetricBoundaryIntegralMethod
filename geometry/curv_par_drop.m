function [k,p,d1,d2,r,dtheta]=curv_par_drop(a,b,n,m,j,q,xc,yc)

  %preallocation
  p=zeros(1,q+1);
  r=zeros(1,q+1);
  dtheta = zeros(1,q);
  %d1=zeros(1,q+1);
  %d2=zeros(1,q+1);
  k=zeros(1,q);
  
  %solo in questo caso
  %xc = 2;
  %yc = 0;
  
  %compute radius
  for i=1:q+1
      aaa = a(n+m+j+i+1);
      bbb = b(n+m+j+i+1);
      r(i) = sqrt((aaa-xc)^2+(bbb-yc)^2);
  end
  
  %compute delta theta
  for i=1:q
      a1 = a(n+m+j+i+1);
      a2 = a(n+m+j+i+2);
      dx1 = abs(a1-xc);
      dx2 = abs(a2-xc);
      theta1 = acos(dx1/r(i));
      theta2 = acos(dx2/r(i+1));
      dtheta(i) = abs(theta1-theta2);
  end

  %calculate the parallel component of the interface curvature with central
  %difference in the nodes (first order for second derivative) in non uniform grid
%   for i=2:q
%       dx1=a(n+m+j+i+1)-a(i+n+m+j);
%       dx2=a(n+m+j+i+2)-a(n+m+j+i+1);
%       d1 = (b(i+n+m+j+2)-(1-(dx1/dx2)^2)*b(i+n+m+j+1)-(dx1/dx2)^2*b(i+n+m+j))/(dx1+(dx1^2/dx2));
%       d2 = (b(i+n+m+j+2)-(1+dx1/dx2)*b(i+n+m+j+1)+(dx1/dx2)*b(i+n+m+j))/(dx1/2*(dx1+dx2));
%       p(i)=-d2/(1+d1^2)^1.5;
%   end
  
  for i=2:q
      %dx1=pi/q;
      %dx2=pi/q;
      dx1 = dtheta(i-1);
      dx2 = dtheta(i);
      d1 = (r(i+1)-(1-(dx1/dx2)^2)*r(i)-(dx1/dx2)^2*r(i-1))/(dx1+(dx1^2/dx2));
      d2 = (r(i+1)-(1+dx1/dx2)*r(i)+(dx1/dx2)*r(i-1))/(dx1/2*(dx1+dx2));
      p(i)=abs(r(i)^2+2*d1^2-r(i)*d2^2)/(r(i)^2+d1^2)^1.5;
  end
  
  %use finete difference forward and backward (non uniform grid) for the starting and ending
  %point
  dx1=dtheta(1);
  dx2=dtheta(1)+dtheta(2);
  d1f = (r(2)-(1-(dx1/dx2)^2)*r(1)-(dx1/dx2)^2*r(3))/(dx1*(1-dx1/dx2));
  d2f = (-(dx1/dx2)*r(3)+r(2)-(1-(dx1/dx2))*r(1))/((dx1/2)*(dx1-dx2));
  p(1)=abs(r(i)^2+2*d1f^2-r(i)*d2f^2)/(r(i)^2+d1f^2)^1.5;

  dx1=dtheta(end)+dtheta(end-1);
  dx2=dtheta(end);
  d1b = (r(end-2)-(dx1/dx2)^2*r(end-1)-(1-(dx1/dx2)^2)*r(end))/(dx1*(dx1/dx2-1));
  d2b = (r(end-2)-(dx1/dx2)*r(end-1)-(1-(dx1/dx2))*r(end))/((dx1/2)*(dx1-dx2));
  %p(q+1)=-(-(1-(dx1/dx2))*(b(n+m+j+q+2)-(dx1/dx2)*b(n+m+j+q+1)+b(n+m+j+q))/((dx1/2)*(dx1-dx2)))/(1+((-(1-(dx1/dx2)^2)*b(q+n+m+j+2)-(dx1/dx2)^2*b(q+n+m+j+1)+b(n+m+j+q))/(dx1*(dx1/dx2-1)))^2)^1.5;
  p(q+1)=abs(r(i)^2+2*d1b^2-r(i)*d2b^2)/(r(i)^2+d1b^2)^1.5;
  
  %calculate the curvature in the sing coordinates as average of the nodes
  %curvature
  for i=1:q
      k(i)=(p(i+1)+p(i))/2;
  end
end