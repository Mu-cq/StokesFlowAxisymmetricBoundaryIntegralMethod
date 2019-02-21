%compute Repulsive Force

function F = computeRepulsiveForce(x1,y1,x2,y2,PARAM)

%compute smallest separation
[dist,sign] = geometryRepulsive(x1,y1,x2,y2);
%[dist,~,~,sign] = panelDistance(x,y,panel1,panel2,PARAM);

%critical distance
delta = PARAM.repulsiveOn;
r = dist;

%compute only when surfaces are close
activate = r<=delta;

%compute intensity repulsive force
F = PARAM.coeffRepulsive*(exp(delta-r)-1)/(exp(delta)-1)*activate*sign;