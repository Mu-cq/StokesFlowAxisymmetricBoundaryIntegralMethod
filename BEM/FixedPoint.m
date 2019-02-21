%fixed point function for newton method

function G = FixedPoint(fNonlinear,r,dh)

    J = JacobianHandle(fNonlinear,r,dh);
    u = fNonlinear(r);
    G = r - J\u;

end