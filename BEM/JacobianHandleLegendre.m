%compute jacobian form function handle

function J = JacobianHandleLegendre(f,r0,theta,dh)

    %dh = 1e-05;
    u = f(r0);
    J = zeros(numel(u),numel(r0));
    
    for i = 1:numel(r0)
        
        %series terms
        n = i-1;
        
        %compute legemdre polynomial for perturbing
        %legendre function
        PPP = legendre(n,cos(theta));
        %legendre polynomia for arbitrary shape
        P = PPP(1,:)./(r0'.^2);
        %P = PPP(1,:);
        
        r = r0+dh*P';
        uNow = f(r);
        
        J(:,i) = (uNow-u)./(dh*P');
        
        figure(1)
        plot(r0.*cos(theta'),r0.*sin(theta'))
        grid on
        hold on
        axis equal
        plot(r.*cos(theta'),r.*sin(theta'))
        hold off
        
    end

end