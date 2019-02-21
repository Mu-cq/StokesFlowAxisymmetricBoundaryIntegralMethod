%rescale a function to get a desired integral value

function R = rescaleFunction(x,y,f,df,intFinal)

F = f+df;
R = int_axis_spline_symmetric(x,y,F)-intFinal;