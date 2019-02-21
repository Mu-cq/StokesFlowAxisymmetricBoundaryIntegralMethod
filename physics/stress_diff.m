function [df_x,df_y]=stress_diff(df,r,q)

  %preallocation
  df_x=zeros(1,q);
  df_y=zeros(1,q);

  %compute the x and y component of deltaF
  for i=1:q
      df_x(i)=df(i)*r(1,i);
      df_y(i)=df(i)*r(2,i);
  end
end