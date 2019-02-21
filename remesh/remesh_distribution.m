%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function [c, d, ds, int] = remesh_distribution(a,b,distr)
    
    %compute spline coefficients
    [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(a,b);
    
    %moltiply elements for vectorial opration
    t = 0:0.001:1;
    ttt = repmat(t',1,numel(bx));
    bxxx = repmat(bx,numel(t),1);
    cxxx = repmat(cx,numel(t),1);
    dxxx = repmat(dx,numel(t),1);
    byyy = repmat(by,numel(t),1);
    cyyy = repmat(cy,numel(t),1);
    dyyy = repmat(dy,numel(t),1);
    
    %length of every element
    ds = sum(sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2))*(t(2)-t(1));
    s = sum(ds);    %length of the all line
    curv = [0 cumsum(ds)];  %curvilinear coordinate
    
    int = cumsum(distr);
    %UNIspace = interp1(curv(2:end),int,curv(1):curv(end)/(numel(distr)-1):curv(end),'spline'); %uniform spacing
    space = int/int(end)*s;  %spacing for the nodes with new distribution
    space = space/space(end)*s;  %spacing for the nodes with new distribution
     
    %decide for the spacing between the nodes (based on distribution)
    c = interp1(curv,a,[0 space'],'spline');
    d = interp1(curv,b,[0 space'],'spline');

end