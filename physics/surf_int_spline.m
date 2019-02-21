%compute the surface integral of a solid of rotation around the axis
%approximated with a spline and than using trapezi rule

function [A,dA] = surf_int_spline(x,y)

%     a = find(x==0);
%     b = find(y==0);
%     
%     x = x(1:a(1)-1);
%     y = y(1:b(2));

%     xcm = center_mass(x,y);
%     x = x-xcm;

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
    
%     figure
%     plot(eta,beta)
%     axis equal

    %using trapezi rules for integration
    DX = eta(1:end-1)-eta(2:end);
    DY = beta(1:end-1)-beta(2:end);
    ds = sqrt(DX.*DX+DY.*DY);
    dA = pi*(beta(1:end-1)+beta(2:end)).*ds;
    A = sum(dA);

end