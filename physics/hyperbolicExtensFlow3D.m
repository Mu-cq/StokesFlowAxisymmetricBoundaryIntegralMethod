%extensional flow

function [u,v,w] = hyperbolicExtensFlow3D(x,z,Ca)

    G = Ca/2;
    
    u = -2*G*x;
    v = 0;
    w = 2*G*z;

end