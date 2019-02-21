%compute surface integral with trapezi rule in cylindrical coordiante

function [A,dA] = surf_lailai(x,y)

    %using trapezi rules for integration
    %dtheta = pi/(numel(x)-1);
    
%     a = find(x==0);
%     b = find(y==0);
%     
%     x = x(1:a(1)-1);
%     y = y(1:b(2));
    
    xcm = center_mass(x,y);
    x = x-xcm;
    
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x, y);
    
    M = 10;
    I = numel(x)-1;
    
    eta = zeros(M*I+1,1);
    beta = zeros(M*I+1,1);
    
    t = 0:1/M:0.9;
    
    for i = 1:I
        for k = 1:M
            eta(k+M*(i-1)) = ax(i)+bx(i)*t(k)+cx(i)*t(k)^2+dx(i)*t(k)^3;
            beta(k+M*(i-1)) = ay(i)+by(i)*t(k)+cy(i)*t(k)^2+dy(i)*t(k)^3;
        end
    end
    
    eta(end) = ax(end)+bx(end)+cx(end)+dx(end);
    beta(end) = ay(end)+by(end)+cy(end)+dy(end);
    
    R = sqrt(eta.*eta+beta.*beta);
    
%     figure
%     plot(x,y,'o-')
%     hold on
%     axis equal
%     plot(eta,beta,'or')
%     hold off
    
    theta = acos(eta./R);
    dtheta = theta(2:end)-theta(1:end-1);
    
    sintheta = beta./R;
    dRdtheta = [diff(R)./diff(theta); diff(R(end-1:end))/diff(theta(end-1:end))];
    
    tot = R.*sintheta.*sqrt(dRdtheta.*dRdtheta+R.*R);
    dA = pi*(tot(1:end-1)+tot(2:end)).*dtheta;
    
    A = sum(dA);
    %A = trapz(theta,tot);
    
%     int = zeros(numel(dA),1);
%     for k=1:numel(dA)
%         int(k) = sum(dA(1:k));
%     end
    
%     figure(2)
%     hold on
%     %plot((1:numel(dA))/numel(dA),dA,'o-')
%     plot(theta(2:end),dA,'o-')
%     title('dA')
%     hold off
%     
%     figure(3)
%     hold on
%     plot(theta)
%     title('theta')
%     hold off
%     
%     figure(4)
%     hold on
%     plot(R)
%     title('R')
%     hold off
%     
%     figure(5)
%     hold on
%     plot(eta,beta,'o')
%     axis equal
%     title('domain')
%     hold off
%     
%     figure(6)
%     hold on
%     plot(eta)
%     title('eta')
%     hold off
%     
%     figure(7)
%     hold on
%     plot(beta)
%     title('beta')
%     hold off
%     
%     figure(8)
%     hold on
%     plot(theta,tot,'o')
%     axis equal
%     title('quantity')
%     hold off
%     
%     figure(9)
%     hold on
%     plot(int,'o')
%     title('integral')
%     hold off
    
    

end