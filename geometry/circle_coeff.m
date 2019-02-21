%compute the coefficient of the equation for a circoference given three
%points

function [R,a,b,c] = circle_coeff(x, y)

    x1 = x(1:end-2);
    y1 = y(1:end-2);
    x2 = x(2:end-1);
    y2 = y(2:end-1);
    x3 = x(3:end);
    y3 = y(3:end);

    %known analitically
    b = ((x2.*x2+y2.*y2-x3.*x3-y3.*y3).*(x2-x1)-(x1.*x1+y1.*y1-x2.*x2-y2.*y2).*(x3-x2))./((y1-y2).*(x3-x2)-(y2-y3).*(x2-x1));
    a = (x2.*x2+y2.*y2-x3.*x3-y3.*y3+b.*(y2-y3))./(x3-x2);
    c = -(x3.*x3+y3.*y3+a.*x3+b.*y3);
    
    R = sqrt(-c+a./2.*a./2+b./2.*b./2);
    
%     figure
%     theta = 0:2*pi/100000:2*pi;
%     plot(R(1)*cos(theta)-a(1)/2,R(1)*sin(theta)-b(1)/2,'-')
%     hold on
%     axis equal
%     plot(R(2)*cos(theta)-a(2)/2,R(2)*sin(theta)-b(2)/2,'-r')
%     plot((R(1)+R(2))/2*cos(theta)-(a(1)*R(1)+a(2)*R(2))/((R(1)+R(2))/2)/2,(R(1)+R(2))/2*sin(theta)-(b(1)*R(1)+b(2)*R(2))/((R(1)+R(2))/2)/2,'-k')

end