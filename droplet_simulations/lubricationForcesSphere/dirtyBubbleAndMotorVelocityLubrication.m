%velocity of dirty bubble with lubrication

function [Ub,Um] = dirtyBubbleAndMotorVelocityLubrication(a,ap,theta,gap,L,H)

%compute constants
A = -12*a^1.5*pi*sin(theta)^2/(2*gap)^1.5;
B = 12*a^1.5*pi*sin(theta)*ap/(2*gap)^1.5;
C = -12*pi*sqrt(a)*ap*sin(theta)/(2*gap)^0.5;
D = -4*pi*cos(theta)^2/(2*gap)^0.5;
E = -4*pi*L/(log(L/H)+log(2)-0.5)/(2*pi*(gap+a)*cos(theta));

%new constant
C1 = (B+C)/(A+D);
C2 = (B+C)/(A+D+E);
C3 = (A+D)/(A+D+E);

%bubble velocity
Ub = (C2-C1)/(1-C3);

%error('Result is not correct')

%motor velocity
Um = (C2-C1)/(1-C3)+C1;