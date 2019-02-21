

function tParametricHere = splitElements(tParametricBase,panel1,panel2,PARAM)

tParametricHere = tParametricBase{panel1};
tParametric = tParametricBase;

%get necessary refinement
maxLevelRefine = ceil(log(PARAM.maxElem(panel1)/PARAM.minSizeElemRemesh(panel1))/log(2));

for k = 1:maxLevelRefine
    
    tParametric{panel1} = tParametricHere;
    
    %build new shape
    [xCell,yCell] = buildGeometryPanelsParametric(tParametric,PARAM);
    x = xCell{panel1};
    y = yCell{panel1};
    xOther = xCell{panel2};
    yOther = yCell{panel2};

    count = 1;
    for i = 1:numel(x)-1
        
       %element coordinates
       x1 = x(i);
       x2 = x(i+1);
       y1 = y(i);
       y2 = y(i+1);
       xm = (x1+x2)/2;
       ym = (y1+y2)/2;
       dl = sqrt((x1-x2)^2+(y1-y2)^2);

       %check from middle point to other shape
       if PARAM.geometryPanel(panel2)==0  % is straight line

            if xm>min(xOther) && xm<max(xOther)
                
            %compute minimum distance between line and point
            mLine = (yOther(end)-yOther(1))/(xOther(end)-xOther(1));
            aLine = mLine;  bLine = -1; cLine = yOther(1)-xOther(1)*mLine;
            distHere = abs(aLine*xm+bLine*ym+cLine)/sqrt(aLine^2+bLine^2);
            
            else
            
            distHere = min(sqrt((xm-xOther).^2+(ym-yOther).^2));
            
            end
           
       elseif PARAM.geometryPanel(panel2)==1  % is an arc
               
            %compute minimum distance between sphere and point
            xcmBubble = PARAM.x0_Circle(panel2);
            ycmBubble = PARAM.y0_Circle(panel2);
            rBubble = PARAM.rArc(panel2);
            distHere = sqrt((xcmBubble-xm)^2+(ycmBubble-ym)^2)-rBubble;
            
       elseif PARAM.geometryPanel(panel2)==2  % is arbitrary shape, like drop or bubble
           
            ppx = PARAM.ppx(panel2);
            ppy = PARAM.ppy(panel2);
            tHere = linspace(0,1,100);
            xHere = ppval(ppx,tHere);
            yHere = ppval(ppy,tHere);
               
            %compute minimum distance between sphere and point
            distHere = sqrt((xHere-xm).^2+(yHere-ym).^2);
            distHere = min(distHere);
               
       end

       if dl>distHere/PARAM.coeffDist

           if PARAM.geometryPanel(panel1)==0 % is a straight line

                %new node
                newNodeParametric = (tParametricHere(count+1)+tParametricHere(count))/2;
                tParametricHere = sort([tParametricHere newNodeParametric]);
            
            elseif PARAM.geometryPanel(panel1)==1 % is an arc

                %new node
                newNodeParametric = (tParametricHere(count+1)+tParametricHere(count))/2;
                tParametricHere = sort([tParametricHere newNodeParametric]);

            else

                error('Not implemented')

            end
        
        %new node position
        count = count+1;

       end
       
       count = count+1;

    end

end

