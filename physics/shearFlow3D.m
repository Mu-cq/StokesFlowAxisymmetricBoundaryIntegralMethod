%extensional flow

function [u,v,w] = shearFlow3D(x,Ca)

    G = Ca/2;
    
    u = G*x;
    v = 0;
    w = 0;

end