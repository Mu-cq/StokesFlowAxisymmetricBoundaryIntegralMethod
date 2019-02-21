%extensional flow

function [u,v] = extens_flow(x,y,Ca)

    G = Ca/2;
    
    u = 2*G*x';
    v = -G*y';

end