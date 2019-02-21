%compute slip velocity once he species concentration are known

function l = computeArcLenghtBlock(x,y,PARAM,block)

%change variables name
PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

%get block coordinates
[x,y,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,block);

if PARAM.orderGeometry(panelRange(1))==0
    
    if numel(PARAM.orderGeometry(panelRange))>1
        if PARAM.orderGeometry(panelRange(1))~=PARAM.orderGeometry(panelRange(2))
            error('Panels on the same block need same geometrical order')
        end
    end
       
    %compute arc lenght constant elem
    dl = sqrt(diff(x).^2+diff(y).^2);
    l = zeros(numel(dl)-1,1);
    l(1) = dl(1)/2;
    for i = 2:numel(x)-1
        
        l(i) = (dl(i-1)+dl(i))/2 + l(i-1);
        
    end
       
elseif PARAM.orderVariable(panelRange(1))==0 && PARAM.orderGeometry(panelRange(1))==1     % in the midlle of the curved element
    
    error('Not implemented')
    
elseif PARAM.orderVariable(panelRange(1))==1 && (PARAM.orderGeometry(panelRange(1))==0 || PARAM.orderGeometry(panelRange(1))==1)     % in the node (with splines)
    
    %compute spline coeff
    if PARAM.SPlinesType(panelRange(1))==1
                    [ax,bx,cx,dx,ay,by,cy,dy] = spline_natural(x,y);
    elseif PARAM.SPlinesType(panelRange(1))==2
                    [ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(x,y);
    end
                
    t = 0:0.05:1;
    ttt = repmat(t',1,numel(ax));
    axxx = repmat(ax,numel(t),1);
    bxxx = repmat(bx,numel(t),1);
    cxxx = repmat(cx,numel(t),1);
    dxxx = repmat(dx,numel(t),1);
    ayyy = repmat(ay,numel(t),1);
    byyy = repmat(by,numel(t),1);
    cyyy = repmat(cy,numel(t),1);
    dyyy = repmat(dy,numel(t),1);
    
    xSP = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    ySP = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    dl = sum(sqrt(diff(xSP).^2+diff(ySP).^2));
    l = [0 cumsum(dl)];
       
end









