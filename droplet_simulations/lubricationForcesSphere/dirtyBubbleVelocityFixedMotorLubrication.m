%velocity of dirty bubble with lubrication

function Ub = dirtyBubbleVelocityFixedMotorLubrication(a,ap,theta,gap)

%compute constants
A = -12*a^1.5*pi*sin(theta)^2/(2*gap)^1.5;
B = 12*a^1.5*pi*sin(theta)*ap/(2*gap)^1.5;
C = -12*pi*sqrt(a)*ap*sin(theta)/(2*gap)^0.5;
D = -4*pi*cos(theta)^2/(2*gap)^0.5;

%bubble velocity
Ub = -(B+C)/(A+D);