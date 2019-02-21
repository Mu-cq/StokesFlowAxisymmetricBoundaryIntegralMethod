%compute jacobian form function handle

function [J,DH] = JacobianHandleAlpha(f,r0,dh)

    %dh = 1e-05;
    u = f(r0);
    J = zeros(numel(u),numel(r0));
    DH = zeros(numel(u),numel(r0));
    
    for i = 1:numel(r0)
        
        r = r0; r(i) = r0(i)+dh;
        uNow = f(r);
        
        %dhHere = rNew-r0;
        %DH(:,i) = dhHere;
        
        %J(:,i) = (uNow-u)./dhHere;
        J(:,i) = (uNow-u)/dh;
        %J = J + (uNow-u)/dhHere;
        
    end

end