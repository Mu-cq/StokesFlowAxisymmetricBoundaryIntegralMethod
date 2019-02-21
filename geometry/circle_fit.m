%compute the circonference coefficient given two points and the radius (in fact I will find two set of coefficient)

function [a1, b1, c1, a2, b2, c2] = circle_fit(x, y, r)

    x1 = x(1:end-1);
    y1 = y(1:end-1);
    x2 = x(2:end);
    y2 = y(2:end);
    
    %distance between the two points
    ds = sqrt((x1-x2).*(x1-x2)+(y1-y2)*(y1-y2));
    
    if sum(ds>2*r)
        display('Warning: at least one radius is too small to reach the desired points')
    end

    f1 = 0.5*(x1.*x2.*(x1-x2)+y1.*y1.*x2-y2.*y2.*x1)./(y2.*x1-y1.*x2);
    f2 = 0.5*(x2-x1)./(y2.*x1-y1.*x2);
    %f3 = 0.5*(x1.*x1+y1.*y1)./x1+0.5*(x1.*x2.*(x1-x2)+y1.*y1.*x2-y2.*y2.*x1)./(x1.*(y2.*x1-y1.*x2));
    f3 = 0.5*(x1.*x1+y1.*y1+2*f1.*y1)./x1;
    %f4 = 0.5*(x2-x1)./(x1.*(y2.*x1-y1.*x2));
    f4 = 0.5*(2*f2.*y1+1)./x1;
    
    alpha = f4.*f4+f2.*f2;
    beta = 2*f3.*f4+2*f1.*f2-1;
    gamma = f1.*f1+f3.*f3-r.*r;
    
    c1 = (-beta+sqrt(beta.*beta-4*alpha.*gamma))./(2*alpha);
    a1 = -2*f3-2*c1.*f4;
    b1 = 2*f1+2*c1.*f2;
    
    c2 = (-beta-sqrt(beta.*beta-4*alpha.*gamma))./(2*alpha);
    a2 = -2*f3-2*c2.*f4;
    b2 = 2*f1+2*c2.*f2;
    
    
    figure
    theta = 0:2*pi/100:2*pi;
    plot(r*cos(theta)-a1/2,r*sin(theta)-b1/2,'o-')
    hold on
    axis equal
    plot(r*cos(theta)-a2/2,r*sin(theta)-b2/2,'o-r')

end