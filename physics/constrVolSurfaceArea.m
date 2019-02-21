%build new shape from stretching a base shape and enforging geometrical constrint (here max lenght and volume)

function [DisEquality,Equality] = constrVolSurfaceArea(x,y,alpha,Atarget,V0)
    
    %stretch shape
    x = alpha(1)*x;
    y = alpha(2)*y;
    
    %residual
    DisEquality = [];
    A = surf_gauss_vect(x,y);
    V = axis_int_gauss_vect(x,y);
    Equality = [A-Atarget V-V0];

end