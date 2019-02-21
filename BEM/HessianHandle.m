%compute jacobian form function handle

function H = HessianHandle(f,r0)

    dh = 1e-06;
    u = f(r0);
    H = zeros(numel(r0),numel(r0),numel(u));
    
    for i = 1:numel(r0)
        
        r = r0; r(i) = r0(i)+dh;
        uNow1 = f(r);
        
        for k = 1:numel(r0)
        
%             r = r0; r(k) = r0(k)+dh; r(i) = r(i) + dh;
%             uNow2 = f(r);
%         
%             H(i,k,:) = (uNow2-2*uNow1+u)/dh^2;
            
            r = r0; r(i) = r0(i)+dh; r(k) = r0(k) - dh;
            uNow2 = f(r);
        
            H(i,k,:) = (uNow2-2*u+uNow1)/dh^2;
        
        end
        
    end

end