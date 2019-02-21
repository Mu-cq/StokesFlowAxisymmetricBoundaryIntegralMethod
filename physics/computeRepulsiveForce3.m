%compute Repulsive Force

function F = computeRepulsiveForce3(x1,y1,x2,y2,PARAM)

%compute smallest separation
[dist,sign] = geometryRepulsive(x1,y1,x2,y2);

%critical distance
delta = PARAM.debeyLenght;
r = dist;

%compute intensity repulsive force
F = PARAM.coeffRepulsive*exp(-delta*r)*sign;