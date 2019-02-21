%build new shape from stretching a base shape and enforging geometrical constrint (here max lenght and volume)

function R = extrapolateOnlyLength(x,y,alpha,fTarget)
    
    %stretch shape
    x = alpha(1)*x;
    y = alpha(2)*y;
    
    %residual
    xcm = center_mass(x,y);
    f = norm(sqrt((x-xcm).^2+y.^2)-1,inf);
    R(1) = abs(f-fTarget);

end