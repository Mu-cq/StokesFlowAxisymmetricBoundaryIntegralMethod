%disjoining pressure

function PI = disj_pressure(h,stop,stretchY)

%normalize
h = h/stop;

%LDVO

%van der Walls
C1 = 5*10^-2;
PIm = C1./h.^2;

%electrostaic
C2 = -50;
C3 = 7;
PIe = C2*exp(-C3*h);

%total value of disjoining pressure
PI = PIm+PIe;

%modify intensity and shut down when h>stop
PI = stretchY*PI.*(h<=1);

end

