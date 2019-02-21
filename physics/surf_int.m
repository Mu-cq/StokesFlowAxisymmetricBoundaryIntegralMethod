%compute the surface integral of a solid of rotation around the axis

function A = surf_int(x,y)

    %using trapezi rules for integration
    dx = (x(1:end-1)-x(2:end));
    dy = (y(1:end-1)-y(2:end));
    dA = sqrt(dx.*dx+dy.*dy).*(y(1:end-1)+y(2:end));
    A = pi*sum(dA);

end