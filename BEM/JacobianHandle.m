%compute jacobian form function handle

function J = JacobianHandle(f,r0,dh)

    u = f(r0);
    J = zeros(numel(u),numel(r0));
       
    for i = 1:numel(r0)
        
        r = r0; r(i) = r0(i)+dh;
        uNow = f(r);
        
        J(:,i) = (uNow-u)/dh;
                
    end

end