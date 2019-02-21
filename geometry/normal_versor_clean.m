function [N,r]=normal_versor_clean(a,b,q)

  %obtains the normal of the inteface elements (pointing upward)
  r=zeros(2,q);
  N=zeros(2,q+1);
  
  %versor normal to the singularity
  for i=1:q;
      alfa=(b(i)-b(i+1))/(a(i)-a(i+1));
      norma=sqrt(alfa^2+1);
      
      %problem in case of strong bending of the interface, use vector
      %product to figureout the signs of the normal vectors
      n=[-alfa/norma 1/norma 0]';
      %tangent vector
      t = [a(i+1)-a(i) b(i+1)-b(i) 0]';
      check = cross(n,t);
      
      if check(3)>=0
        r(:,i)=[-alfa/norma 1/norma]';
      else
          r(:,i)=[alfa/norma -1/norma]';
      end
  end
  
  %obtain the normal of the nodes of the interface as average of the normal
  %of the 2 adjacent elements
  for i=2:q
      N(:,i)=(r(:,i)+r(:,i-1))/2;
  end
  
  N(:,1) = r(:,1);
  N(:,q+1) = r(:,q);
  
end