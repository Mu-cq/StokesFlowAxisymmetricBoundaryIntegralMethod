%compute weights for SPlines perpendicula to the axis

function [weightA,weightB,weightC,weightD] = splinesWeightsSymmetric(n)

    error('Symmetric spline seems not perfect close to rhe axis')

    %matrix which sets the kind of splines
    A = diag(4*ones(1,n));
    A(1,1) = 2;
    A(end,end) = 2;
    A = A+diag(1*ones(1,n-1),1);
    A = A+diag(1*ones(1,n-1),-1);
    
    %if I want symmetry condition
    A(1,1) = 1;
    A(1,2) = 0;
    A(end,end-1) = 0;
    A(end,end) = 1;

    %a coefficient
    weightA = diag(ones(n,1));
    weightA = weightA(1:end-1,:);
    
    %b coefficient
    B = diag(1*ones(1,n-1),1) - diag(1*ones(1,n-1),-1);
    B(1,1) = -1;    B(end,end) = 1;
    cut = 3*(A\B);
    weightB = cut(1:end-1,:);
    
    %c coefficient
    D = A\B;
    Di = 3*D(1:end-1,:); Di1 = 3*D(2:end,:);
    B = -diag(ones(1,n)) + diag(ones(1,n-1),1);
    B = B(1:end-1,:);
    weightC = 3*B -2*Di -Di1;
    
    %d coefficient
    B = diag(ones(1,n)) - diag(ones(1,n-1),1);
    B = B(1:end-1,:);
    weightD = 2*B + Di + Di1;

end