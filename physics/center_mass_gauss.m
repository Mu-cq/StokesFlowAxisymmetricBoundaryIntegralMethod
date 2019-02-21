%compute center of mass of a solid of rotation

function xcm = center_mass_gauss(x,y)

error('bug')

        %splines coefficients
        [~, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(x,y);

        %compute normal vector
        N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
            -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

        %compute numerator using divergence theorem
        fNum = sum([x.^2/2; zeros(1,numel(x))].*N);
        num = int_axis_spline_symmetric_coeff(bx,cx,dx,ay,by,cy,dy,fNum);
        
        %compute denominator using divergence theorem
        fDen = sum([x/2; zeros(1,numel(x))].*N);
        den = int_axis_spline_symmetric_coeff(bx,cx,dx,ay,by,cy,dy,fDen);
        
        %using divergence theorem
        xcm = num/den;

return