%distance between a series of nodes and a droplet

function dist = NodeDropDist(ax,bx,cx,dx,ay,by,cy,dy,a,b)
    
    %compute splines coordinates
    t = 0:0.1:0.9;
    ttt = repmat(t,1,numel(ax));
    axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
    bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
    cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
    dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
    ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
    byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
    cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
    dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));
        
    %splines coordinates
    x = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    y = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    
    %compute the angular coefficient of the line perpendicular to the spline
    no = sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2);
    nx = (byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2)./no;
    ny = -(bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2)./no;
    mINV = nx./ny;
    
    %here I fix the places where the nearly singular tretament is performed
    %copy for vectorial operation
    xxx = repmat(x',1,numel(a));
    yyy = repmat(y',1,numel(a));
    aaa = repmat(a,numel(x),1);
    bbb = repmat(b,numel(x),1);
    mmm = repmat(mINV',1,numel(a));
    %find the line paerdendicular to the spline curve and passing for the
    %node
    RES = mmm.*(bbb-yyy)-(aaa-xxx);
    [~,ind] = min(abs(RES));
    %ind(2:end) = ind(2:end) + cumsum(ones(1,numel(ind)-1));
    %dist = min(abs(RES));
    
    dist = sqrt((a-xxx(ind)).^2 + (b-yyy(ind)).^2);
    %aMIN = a;
    %aMIN = aaa(ind);
    
    
    
  
end