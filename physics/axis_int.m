%compute the volume integral of a solid of rotation around the axis

function V = axis_int(x,y)

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

    %using trapezi rules for integration
    DX = (eta(1:end-1)-eta(2:end));
    dV = (beta(1:end-1).*beta(1:end-1)+beta(2:end).*beta(2:end)).*DX/2;
    V = pi*sum(dV);

end