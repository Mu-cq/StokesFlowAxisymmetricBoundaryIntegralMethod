%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y] = draw_circle(x1,y1,R,N,wave,ampli)

%create the empty matrix
theta = 0:pi/N:pi;

%perturbed radius
R = R*(1+ampli*cos(wave*theta));

%build the circle (half)
X = R.*cos(theta)+x1;
Y = R.*sin(theta)+y1;

% figure
% plot(R)
% 
% figure
% plot(X,Y,'o-')
end