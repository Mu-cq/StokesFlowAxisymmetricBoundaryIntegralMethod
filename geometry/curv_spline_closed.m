%compute curvature using splines for a closed geometry

function curv = curv_spline_closed(bx,by,cx,cy)

    num = 2*bx.*cy-2*by.*cx;
    den = (bx.*bx+by.*by).^1.5;
    
    %last = ((bx(end)+2*cx(end)+3*dx(end))*(2*cy(end)+6*dy(end))-(by(end)+2*cy(end)+3*dy(end))*(2*cx(end)+6*dx(end)))/...
     %   ((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))^1.5;
    
    %curv = [num./den last];
    curv = num./den;
        
%     figure
%     plot(curv,'o-')

end