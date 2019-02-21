%newton method

function x = myNewtonMethod(fNonlinear,xIN,errTol)

dh = 1e-5;

%newton method loop
count = 1;
R = 1;
x = xIN;
while norm(R,Inf)>errTol
    
    %compute residuals
    R = fNonlinear(x);
    
    %compute Jacobian
    J = JacobianHandle(fNonlinear,x,dh);
    
    %solve linearized system
    dx = -J\R;
    
    %update solution
    x = x+dx;
    
    if count>100
        warning('Newton Method not converging')
        break;
    end
    
    count = count+1;
    
end