%compute the coefficient of the equation for a circoference given three
%points

function [R,a,b,c] = CurvDrop2(x, y, n)

    x1 = x(1:end-2);
    y1 = y(1:end-2);
    x2 = x(2:end-1);
    y2 = y(2:end-1);
    x3 = x(3:end);
    y3 = y(3:end);

    %known analitically
    b = ((x2.*x2+y2.*y2-x3.*x3-y3.*y3).*(x2-x1)-(x1.*x1+y1.*y1-x2.*x2-y2.*y2).*(x3-x2))./((y1-y2).*(x3-x2)-(y2-y3).*(x2-x1));
    a = (x2.*x2+y2.*y2-x3.*x3-y3.*y3+b.*(y2-y3))./(x3-x2);
    c = -(x3.*x3+y3.*y3+a.*x3+b.*y3);
    
    R = sqrt(-c+a./2.*a./2+b./2.*b./2);
    
    %figure out the sign of the curvature
    %radius vector
    r = [x2+a/2; y2+b/2];
    
    
    %scalar product between the normal at the nodes and the radius vector
    
    for i = 1:length(a)
        discr = r(:,i)'*n(:,i);
        %display(discr)

        if discr<0
            R(i) = -R(i);
        end
    
    end
    

end