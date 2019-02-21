%print to screen information when starting the simulation

function printToScreenStokes2d(PARAM)

if numel(PARAM.n)~=sum(PARAM.panels)
    error('Total panels do not correspond to panels in block')
end

%LAPLACE SOLVER
disp('STOKES SOLVER')
if PARAM.STstokes==1
    disp('Singularity treatment for Stokes equation is activated')
elseif PARAM.STstokes==0
    disp('Singularity treatment for Stokes equation is NOT activated')
end

for i = 1:numel(PARAM.panels)
    
    if i==1
        sumPanels = 0;
    else
        sumPanels = sum(PARAM.panels(1:i-1));
    end
    for k = 1:PARAM.panels(i)
       
        numPan = k+sumPanels;
        if PARAM.typeBCstokes(numPan)==0
            BC = 'prescribed horizontal and vertical velocity';
        elseif PARAM.typeBCstokes(numPan)==1
            BC = 'prescribed normal velocity (tangent is zero)';
        elseif PARAM.typeBCstokes(numPan)==2
            BC = 'prescribed normal stress (tangent is zero)';
        elseif PARAM.typeBCstokes(numPan)==3
            BC = 'prescribed tangent velocity';
        elseif PARAM.typeBCstokes(numPan)==4
            BC = 'prescribed horizontal velocity';
        elseif PARAM.typeBCstokes(numPan)==5
            BC = 'prescribed normal velocity and tangent stress';
        elseif PARAM.typeBCstokes(numPan)==6
            BC = 'rigid body motion';
        elseif PARAM.typeBCstokes(numPan)==7
            BC = 'prescribed horizontal stress and vertical velocity';
        elseif PARAM.typeBCstokes(numPan)==8
            BC = 'prescribed normal and tangent stress';
        elseif PARAM.typeBCstokes(numPan)==9
            BC = 'prescribed normal stress and tangent velocity';
        else
            error('BC not implemented')
        end
        
        if PARAM.panelType(numPan)==0
            wallType = 'fixed wall';
        elseif PARAM.panelType(numPan)==1
            wallType = 'moving wall';
        elseif PARAM.panelType(numPan)==2
            wallType = 'droplet interface';
        end
        
        if PARAM.orderGeometryStokes(numPan)==0
            orderGeo = 'straight';
        elseif PARAM.orderGeometryStokes(numPan)==1
            orderGeo = 'curved';
        end
        
        if PARAM.orderVariableStokes(numPan)==0
            orderVar = 'P0';
        elseif PARAM.orderVariableStokes(numPan)==1
            orderVar = 'P1';
        end
        
        disp(['Panel ' num2str(k) ' of block ' num2str(i) ' is initialized as a ' wallType ' with ' num2str(PARAM.n(numPan)) ' ' orderVar ' ' orderGeo ' elements, BC is ' BC])
       
    end
    
end

count = 0;
for i = 1:numel(PARAM.panels)
    
   for k = 1:PARAM.panels(i)
       
      if PARAM.blockType(i)~=PARAM.panelType(count+k)
          error('Panel has to be of the same type as the block they belong to')
      end
       
   end
   count = count+k;
    
end

if PARAM.addFlow==1
    display(['Underlying velocity is U=' num2str(PARAM.Uunder)])
elseif PARAM.addFlow==2
    display(['Underlying velocity is an extensional flow Ca=' num2str(PARAM.Ca)])
end

for i = 1:numel(PARAM.panels)
    if PARAM.deflationBlock(i)==1
        display(['Deflation is applied on block ' num2str(i)])
    else
        display(['NO deflation is applied on block ' num2str(i)])
    end
end

if sum(PARAM.repulsiveForces)>0
    
   disp('Repulsive forces are active') 
   
else
    
    disp('Repulsive forces are NOT active')
    
end






