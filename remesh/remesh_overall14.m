%remesh depending on the droplet volume (the distribtion is kept good with mesh stabilization)

function [c, d] = remesh_overall11(a,b)

    %compute splines coordinates
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
    
    c = zeros(1.1*numel(ax)+1,1);
    d = zeros(1.1*numel(ax)+1,1);
    
    c(1) = a(1);
    d(1) = b(1);
    
    for i = 2:2:2*numel(ax)
        
        t = 0.5;
        c(i) = ax(i/2) + bx(i/2)*t + cx(i/2)*t^2 + dx(i/2)*t^3;
        d(i) = ay(i/2) + by(i/2)*t + cy(i/2)*t^2 + dy(i/2)*t^3;
        
        t = 1;
        c(i+1) = ax(i/2) + bx(i/2)*t + cx(i/2)*t^2 + dx(i/2)*t^3;
        d(i+1) = ay(i/2) + by(i/2)*t + cy(i/2)*t^2 + dy(i/2)*t^3;
        
    end

end