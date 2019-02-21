%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function [c, d, ds, int] = remeshFixedGeometry(a,b,distr,PARAM)
    
    %length of every element
    ds = 2*PARAM.L;
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