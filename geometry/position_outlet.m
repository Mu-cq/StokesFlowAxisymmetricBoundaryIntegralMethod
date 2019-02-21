function I=position_outlet(Ysing,b,n)

  %preallocation
  l=100*ones(1,n);

  %find interface position at the outlet
  for i=1:n
      if Ysing(i)>=b(length(b))
          l(i)=Ysing(i)-b(length(b));
      end
  end
  
  [~,I]=min(l);
  
end