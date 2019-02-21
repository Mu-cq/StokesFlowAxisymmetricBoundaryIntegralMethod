function [N,r]=normal_versor(a,b,n,m,j,q)

  %obtains the normal of the inteface elements (pointing upward)
  r=zeros(2,q);
  N=zeros(2,q);
  
  %versor normal to the singularity
  for i=1:q;
      alfa=(b(i+n+m+j+2)-b(i+n+m+j+1))/(a(i+n+m+j+2)-a(i+n+m+j+1));
      norma=sqrt(alfa^2+1);
      r(:,i)=[-alfa/norma 1/norma]';
  end
  
  %obtain the normal of the nodes of the interface as average of the normal
  %of the 2 adjacent elements
  for i=1:q-1
      N(:,i)=(r(:,i)+r(:,i+1))/2;
  end
  
  N(:,q)=r(:,q);
  
end