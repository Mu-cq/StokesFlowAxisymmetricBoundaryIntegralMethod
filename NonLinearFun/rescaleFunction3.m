%rescale a function to get a desired integral value

function R = rescaleFunction3(x,y,df,Const,intFinal,nx)

F = (df + Const.*abs(nx))';
R = int_axis_spline_symmetric(x,y,F)-intFinal;