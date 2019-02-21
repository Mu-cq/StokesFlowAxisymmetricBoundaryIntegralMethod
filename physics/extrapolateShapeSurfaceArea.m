%build new shape from stretching a base shape and enforging geometrical constrint (here max lenght and volume)

function R = extrapolateShapeSurfaceArea(x,y,alpha,Atarget,V0)
    
    %stretch shape
    x = alpha(1)*x;
    y = alpha(2)*y;
    
    %residual
    A = surf_gauss_vect(x,y);
    R(1) = A-Atarget;
    R(2) = axis_int_gauss_vect(x,y)-V0;

end