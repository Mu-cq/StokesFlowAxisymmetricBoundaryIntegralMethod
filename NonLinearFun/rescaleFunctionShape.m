%rescale a shape to get a desired integral value

function R = rescaleFunctionShape(x,y,gamma,alpha,intFinal)

%rescale x coordinate
y = alpha*y;

%compute stresses in the axial direction
dfX = normalStresses(x,y,gamma);

%compute residuals
R = int_axis_spline_symmetric(x,y,dfX')-intFinal;