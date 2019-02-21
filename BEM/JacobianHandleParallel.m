%compute jacobian form function handle

function J = JacobianHandleParallel(f,r0,dh)

    u = f(r0);
    J = zeros(numel(u),numel(r0));
    
%     if numel(r0)~=numel(u)
%         warning('Non squared Jacobian')
%     end
    
    %parfor
    parfor (i = 1:numel(r0))

        r = r0; r(i) = r0(i)+dh;
        uNow = f(r);
        
        J(:,i) = (uNow-u)/dh;
        
    end

end