%given the coordinate of the end points of the elements, compute the
%coefficient of the arc passing from two nodes before and one after and the
%angles identifing the midpoint of the element

function [theta,dtheta,R,a,b,c] = CompArcElem(x,y)

        %compute the coefficient and radius of a circle passing for three
        %points
        [R,a,b,c] = circle_coeff([x(2) x], [-y(2) y]);
        %[R,a,b,c] = circle_coeff([x x(end-1)], [y -y(end-1)]);
        
        %compute the angle corresponding to each element
        dl = sqrt((x(1:end-1)-x(2:end)).*(x(1:end-1)-x(2:end))+(y(1:end-1)-y(2:end)).*(y(1:end-1)-y(2:end)));
        dtheta = 2*asin(dl/2./R);
        
        %compute the angle identifing the midpoint of each arc element
        Ax = -a/2+R;
        Ay = -b/2;
        l1 = sqrt((x(2:end)-Ax).*(x(2:end)-Ax)+(y(2:end)-Ay).*(y(2:end)-Ay));
        l2 = sqrt((x(1:end-1)-Ax).*(x(1:end-1)-Ax)+(y(1:end-1)-Ay).*(y(1:end-1)-Ay));
        check = zeros(2,numel(x)-1);
        for i = 2:numel(x)
            r1 = [R(i-1); 0; 0];
            r2 = [x(i)+a(i-1)/2; y(i)+b(i-1)/2; 0];
            r3 = [x(i-1)+a(i-1)/2; y(i-1)+b(i-1)/2; 0];
            vect1 = cross(r1,r2);
            vect2 = cross(r1,r3);
            check(1,i-1) = vect1(3);
            check(2,i-1) = vect2(3);
        end
        
%         check1 = l1/2./R;
%         check2 = l2/2./R;
        
        %safenet because it is very delicate when l2~2*R, l2 can becomes
        %greater than 2*R due to numerical errors
        for i = 1:numel(l1)
            if (l1(i)>2*R(i))
                l1(i) = 2*R(i);
            end
            if (l2(i)>2*R(i))
                l2(i) = 2*R(i);
            end
        end
        
        theta1 = (check(1,:)>=0).*(2*asin(l1/2./R))+(check(1,:)<0).*(2*pi-2*asin(l1/2./R));%+(check(1,:)<0.0000000001.*check(1,:)>-0.0000000001).*pi;
        theta2 = (check(2,:)>=0).*(2*asin(l2/2./R))+(check(2,:)<0).*(2*pi-2*asin(l2/2./R));%+(l2>2*R).*(pi);
        
        %theta1 is defined greater thean theta2, if one element cross
        %theta=0 this creates problem, exception
        for i = 1:numel(theta1)
            if (theta1(i)+1.5*pi<theta2(i))
                theta1(i) = theta1(i)+2*pi;
            end
        end
        
        theta = (theta1+theta2)/2;
        %exception for theta(1) because when l1==0 it doesn't know where to
        %take +pi or not
        if (real(theta2(1)) < eps)
            theta(1) = theta1(1)/2;
        end
        
        %theta = real(theta);
                

return