%rescale a function to get a desired integral value

function R = rescaleFunction2(x,y,df,Const,intFinal,nx)

F = (df + Const.*nx.^2)';
R = int_axis_spline_symmetric(x,y,F)-intFinal;