%extensional flow

function [u,v] = extens_flow2(x,y,Ca)

    G = Ca/2;
    
    u = 2*G*x;
    v = -G*y;

end