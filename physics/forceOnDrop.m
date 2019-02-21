%compute force acting on droplet, for a droplet shuold be zero by
%definition

function DropForce = forceOnDrop(a,b,K,gamma)

    %compute normal vector
    [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(a', b');
    N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2);...
        -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2)];
    
    %compute force in axial direction, the one in the radial direction is
    %zero by definition
    dfX = K.*N(1,:)'*gamma;
    
    %integration of force in the axial direction
    DropForce = int_axis_spline_symmetric(a',b',dfX);

end