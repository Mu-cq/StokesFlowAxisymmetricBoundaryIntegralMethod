function L=position_inlet(Ysing,R1,n,m,j)

  %preallocation
  l=10*ones(1,n+m+j);

  %find interface position at the inlet
  for i=n+m+1:n+m+j
      if Ysing(i)<=R1
          l(i)=-Ysing(i)+R1;
      end
  end
  
  [~,L]=min(l);
  
end