

function tParametricHere = remeshClosestElements(numelNodes,tParametricBase,indNode,panel1,PARAM,distance)

tParametricHere = tParametricBase{panel1};

%remesh closest two closest elements
if indNode==1 || indNode==numelNodes
        error('Not implemented')
else
        
    ind1 = indNode-1;
    ind2 = indNode;
    ind3 = indNode+1;
    
    dl1 = 1;
    dl2 = 1;
    count = 0;
    while dl1>distance/2
    
        if PARAM.geometryPanel(panel1)==0 % is a straight line

            %new node
            newNodeParametric = (tParametricHere(ind2)+tParametricHere(ind1))/2;
            newDT = tParametricHere(ind2)-newNodeParametric;
            tParametricHere = sort([tParametricHere newNodeParametric]);
            
            %new dl
            dx = PARAM.xEnd(panel1)-PARAM.xStart(panel1);
            dy = PARAM.yEnd(panel1)-PARAM.yStart(panel1);
            dl1 = sqrt(dx^2+dy^2)*newDT;
            
        elseif PARAM.geometryPanel(panel1)==1 % is and arc

            %new node
            newNodeParametric = (tParametricHere(ind2)+tParametricHere(ind1))/2;
            tParametricHere = sort([tParametricHere newNodeParametric]);
            
            %new dl
            dt = tParametricHere(ind2+1)-tParametricHere(ind1+1);
            dl1 = PARAM.rArc(panel1)*dt;
            
        else
            
            error('Not implemented')

        end
        
        %new node position
        ind3 = ind3+1;
        ind2 = ind2+1;
        ind1 = ind1+1;
        
        if count>100
            error('Too many iterations')
        end
        
        if dl1<PARAM.minSizeElemRemesh(panel1)
            warning('Stop remesh because smallest elment is becoming too small')
            break
        end
        
        count = count+1;
    
    end
    
    %new indeces
    %ind2 = indNode+count+1;
    %ind3 = indNode+count+2;
    
    count = 0;
    while dl2>distance/2
    
        if PARAM.geometryPanel(panel1)==0 % is a straight line

            %new node
            newNodeParametric = (tParametricHere(ind3)+tParametricHere(ind2))/2;
            newDT = newNodeParametric-tParametricHere(ind2);
            tParametricHere = sort([tParametricHere newNodeParametric]);
            
            %new dl
            dx = PARAM.xEnd(panel1)-PARAM.xStart(panel1);
            dy = PARAM.yEnd(panel1)-PARAM.yStart(panel1);
            dl2 = sqrt(dx^2+dy^2)*newDT;

        elseif PARAM.geometryPanel(panel1)==1 % is and arc

            %new node
            newNodeParametric = (tParametricHere(ind3)+tParametricHere(ind2))/2;
            tParametricHere = sort([tParametricHere newNodeParametric]);
            
            %new dl
            dt = tParametricHere(ind3)-tParametricHere(ind2);
            dl2 = PARAM.rArc(panel1)*dt;
            
        else
            
            error('Not implemented')

        end
        
        if count>100
            error('Too many iterations')
        end
        
        if dl2<PARAM.minSizeElemRemesh(panel1)
            warning('Stop remesh because smallest elment is becoming too small')
            break
        end
        
        count = count+1;
    
    end
    
end

