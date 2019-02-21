%newton method

function x = myBisectionMethod(fNonlinear,range,errTol,opt,output)

rangeLow = range(1);
rangeHigh = range(2);

%bisection loop
count = 1;
R = 1;
while norm(R,Inf)>errTol
    
    %guess
    x = (rangeHigh+rangeLow)/2;
    
    %compute residuals
    R = fNonlinear(x);
    
    %bisection
    if opt == 1
        
        if R>0

            rangeHigh = x;

        elseif R<0

            rangeLow = x;

        end
    
    elseif opt == 2
        
        if R<0

            rangeHigh = x;

        elseif R>0

            rangeLow = x;

        end
        
    else
        
        error('Invalid option')
        
    end
    
    if count>1e3
        warning('Bisection Method not converging')
        break;
    end
    
    if output==1
        
       disp(['Bisection iteration ' num2str(count) ': residuals=' num2str(norm(R,Inf))])
        
    end
    
    count = count+1;
    
end