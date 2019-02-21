%build new shape from stretching a base shape and enforging geometrical constrint (here max lenght and volume)

function R = extrapolateShape(x,y,alpha,fTarget,V0)
    
    %stretch shape
    x = alpha(1)*x;
    y = alpha(2)*y;
    
    %residual
    xcm = center_mass(x,y);
    f = norm(sqrt((x-xcm).^2+y.^2)-1,inf);
    R(1) = f-fTarget;
    R(2) = axis_int_gauss_vect(x,y)-V0;

end