%finf distance between two lines

function h = findHeight2(xDown,yDown,xUp,yUp)

h = zeros(1,numel(xDown));
for i = 1:numel(xDown)
    
    %two closets points
    [~,ind1] = min(abs(xDown(i)-xUp));
    
    if ind1==1
        ind2 = 2;
    elseif ind1==numel(xUp);
        ind2 = numel(xUp)-1;
    else
        [~,ind2] = min(abs(xDown(i)-[xUp(1:ind1-1) Inf xUp(ind1+1:end)]));
    end
    %interpolate
    x1 = xUp(ind1); y1 = yUp(ind1);
    x2 = xUp(ind2); y2 = yUp(ind2);
    m = (y1-y2)/(x1-x2);
    
    %h(i) = (yUp(ind)-yDown(i));
    h(i) = y1+m*(xDown(i)-x1);
    
end