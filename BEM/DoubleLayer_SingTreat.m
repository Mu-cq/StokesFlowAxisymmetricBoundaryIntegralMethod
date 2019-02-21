%singular treatment for double layer potential

function [A11,A12,A21,A22] = DoubleLayer_SingTreat(A11,A12,A21,A22,T1,T2,D1,D2,choose)

    if choose==1    %normal points outside the closed surface

        DEsing11 = diag((-sum(T1,2)-4*pi));
        DEsing12 = diag(-sum(D1,2));
        DEsing21 = diag(-sum(T2,2));
        DEsing22 = diag((-sum(D2,2)-4*pi));

        A11 = A11 + DEsing11;
        A12 = A12 + DEsing12;
        A21 = A21 + DEsing21;
        A22 = A22 + DEsing22;
    
    elseif choose==2    %normal points inside the closed surface
        
        DEsing11 = diag((-sum(T1,2)+4*pi));
        DEsing12 = diag(-sum(D1,2));
        DEsing21 = diag(-sum(T2,2));
        DEsing22 = diag((-sum(D2,2)+4*pi));

        A11 = A11 + DEsing11;
        A12 = A12 + DEsing12;
        A21 = A21 + DEsing21;
        A22 = A22 + DEsing22;
        
    end

end