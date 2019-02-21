%compute jacobian form function handle

function J = JacobianHandle2ndOrder(f,r0,dh)

    %dh = 1e-05;
    %u = f(r0);
    J = zeros(numel(r0));
    
    for i = 1:numel(r0)
        
        r1 = r0; r1(i) = r0(i)+dh;
        uNow1 = f(r1);
        
        r2 = r0; r2(i) = r0(i)-dh;
        uNow2 = f(r2);
        
        J(:,i) = (uNow1-uNow2)/dh/2;
        
    end

end