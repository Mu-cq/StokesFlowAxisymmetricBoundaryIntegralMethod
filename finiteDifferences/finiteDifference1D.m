%compute the first derivative of acertian order by finite differences 1D
%n is the nuber of nodes and m the order of accuracy

function D = finiteDifference1D(n,m,der)

    D = zeros(n,n,der);
    dt = 1/(n-1);

    if m(1)==2
        
        %coefficient (Wikipedia)
        D(:,:,1) = diag(0.5*ones(n-1,1),1) + diag(-0.5*ones(n-1,1),-1);
    
        %first point
        D(1,[1 2 3],1) = [-1.5 2 -0.5];
        D(end,[end-2 end-1 end],1) = [0.5 -2 1.5];
        
        %divide by dt
        D(:,:,1) = D(:,:,1)/dt;
    
    end
    
    if m(2)==1
        
        %coefficient (Wikipedia)
        D(:,:,2) = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) + diag(-2*ones(n,1));
    
        %first point
        D(1,[1 2 3 4],2) = [2 -5 4 -1];
        D(end,[end-3 end-2 end-1 end],2) = [-1 4 -5 2];
    
        %divide by dt^2
        D(:,:,2) = D(:,:,2)/dt^2;
        
    end

end