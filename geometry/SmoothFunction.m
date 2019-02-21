%smooth function

function y = SmoothFunction(x0,y0)

    %fit smoothing spline and evaluete data where needed
    f = fit(x0,y0,'smoothingspline');
    y = feval(f,x0);

end