%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y] = draw_closed_circle(x1,y1,R,N,wave,ampli)

%create the empty matrix
theta = 0:2*pi/N:2*pi;

%perturbed radius
R = R*(1+ampli*cos(wave*theta));

%build the circle (half)
X = R.*cos(theta)+x1;
Y = R.*sin(theta)+y1;

end