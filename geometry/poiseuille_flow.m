%compute poiseuille profile given flow rate and channel rate

function u = poiseuille_flow(b,Q,R)

    u = 2*Q/pi/R^4*(R^2-b.^2);

end