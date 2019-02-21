%Choose bubble position and size

function [breakNow,initial,tParametricBase,Tsave,PARAM] = initialConditionNextCycleManyBubbles(tParametricBase,T,Y,Tsave,Tend,SaveHowMany,outFun,PARAM)

%number of current bubble
nBubble = numel(PARAM.n)-4;
Y = Y(end,:);

%new t parametric
tParametricBase(1:4) = tParametricBase(1:4);

%break if there is no new nucleation happens
breakNow = 0;
if T(end)==Tsave(end) && nBubble==1
    disp('Cycle ended without nucleating a second bubble')
    breakNow = 1;
end

%next saving time
Tstart = T(end);
Tsave = linspace(Tstart,Tend+Tstart,SaveHowMany);
    
if breakNow==0
    
    %figure out if there is one more bubble or one less
    [~,newBubble,zNew,rNew,indOut] = outFun(0,Y,0);
    
    if newBubble==1
        
        %initial condition for next iteration
        [bubbleLine,ind] = sort([Y(2:nBubble+1) zNew]);
        bubbleRadii = [Y(nBubble+2:end) rNew];
        initial = [Y(1) bubbleLine bubbleRadii(ind)];
        
        %new parameters for new bubble
        PARAM.x0_Circle(5+nBubble) = nan;
        PARAM.y0_Circle(5+nBubble) = 0;
        PARAM.rArc(5+nBubble) = nan;
        
        %parameters for next iteration
        PARAM.n = [PARAM.n(1:4) 50*ones(1,nBubble+1)];
        PARAM.panels = [4 ones(1,nBubble+1)];
        PARAM.blockType = [1 ones(1,nBubble+1)];
        PARAM.remeshProximity = {[] [] 4+(1:1+nBubble) 4+(1:1+nBubble)};
        for i = 1:nBubble+1
           
            PARAM.remeshProximity{i+4} = [3 4 5:i+3 i+5:nBubble+5];
            
            %new t parametric
            tParametricBase{i+4} = linspace(0,pi,PARAM.n(i+4));
            
        end
        PARAM.remesh = [0 0 ones(1,nBubble+3)];
        
    elseif newBubble==-1
        
        if nBubble==1
            error('Bubble cannot all vanish')
        end
        
        %initial condition for next iteration
        bubbleLine = Y(2:nBubble+1);
        bubbleRadii = Y(nBubble+2:end);
        initial = [Y(1) bubbleLine(1:indOut-1) bubbleLine(indOut+1:end) bubbleRadii(1:indOut-1) bubbleRadii(indOut+1:end)];
        
        %parameters for next iteration
        PARAM.n = [PARAM.n(1:4) 50*ones(1,nBubble-1)];
        PARAM.panels = [4 ones(1,nBubble-1)];
        PARAM.blockType = [1 ones(1,nBubble-1)];
        PARAM.remeshProximity = {[] [] 4+(1:nBubble-1) 4+(1:nBubble-1)};
        for i = 1:nBubble-1
           
            PARAM.remeshProximity{i+4} = [3 4 5:i+3 i+5:nBubble+3];
            
            %new t parametric
            tParametricBase{i+4} = linspace(0,pi,PARAM.n(i+4));
            
        end
        PARAM.remesh = [0 0 ones(1,nBubble+1)];
        tParametricBase = tParametricBase(1:4+nBubble-1);
        
    else
        
       error('Not implemented') 
       
    end
        
end
    
    