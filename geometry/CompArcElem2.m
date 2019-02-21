%given the coordinate of the end points of the elements, compute the
%coefficient of the arc passing from two nodes before and one after and the
%angles identifing the midpoint of the element

function [theta,dtheta,R,a,b,c] = CompArcElem2(x,y)

        %compute the coefficient and radius of a circle passing for three
        %points
        [R,a,b,c] = circle_coeff([x(2) x], [-y(2) y]);
        
        l1 = (x(1:end-1)+a/2)./R;
        l2 = (y(1:end-1)+b/2)./R;
        thecos = acos(l1);
        thesin = asin(l2);
        
        theta1 = ((l1 > 0 & l2 > 0)).*(thecos) + ((l1 > 0 & l2 < 0)).*(2*pi-thesin) + ((l1 < 0 & l2 > 0)).*(thecos) + ((l1 < 0 & l2 < 0)).*(2*(pi-thecos)+thecos) + ...
            (l1 == 0 & l2 > 0).*pi/2 + (l1 == 0 & l2 < 0).*(1.5*pi) + (l1 > 0 & l2 == 0).*0 + (l1 < 0 & l2 == 0).*(pi);
        
        l1 = (x(2:end)+a/2)./R;
        l2 = (y(2:end)+b/2)./R;
        thecos = acos(l1);
        thesin = asin(l2);
        
        theta2 = ((l1 >= 0 & l2 >= 0)).*(thecos+thesin)./2 + ((l1 >= 0 & l2 < 0)).*(thesin) + ((l1 < 0 & l2 >= 0)).*(thecos) + ((l1 < 0 & l2 < 0)).*(2*(pi-thecos)+thecos);
        
        %theta2 is defined greater than theat1!
        theta2 = (theta2+pi < theta1).*(theta2+2*pi) + (theta2+pi >= theta1).*(theta2);
        
        dtheta = theta2 - theta1;
        theta = (theta1+theta2)/2;
                

return