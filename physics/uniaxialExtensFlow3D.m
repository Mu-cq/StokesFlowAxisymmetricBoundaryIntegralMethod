%extensional flow

function [u,v,w] = uniaxialExtensFlow3D(x,y,z,Ca)

    G = Ca/2;
    
    u = -G*x;
    v = -G*y;
    w = 2*G*z;

end