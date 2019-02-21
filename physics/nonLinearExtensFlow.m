%extensional flow

function [u,v] = nonLinearExtensFlow(x,y,Ca,Ca2)
    
    u = Ca*x + Ca*Ca2*x.^3;
    v = -0.5*Ca*y - 1.5*Ca*Ca2*y.*x.^2;

end