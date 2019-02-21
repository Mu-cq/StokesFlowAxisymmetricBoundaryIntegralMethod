%Choose bubble position and size

function [breakNow,initial,tParametricBase,Tsave,PARAM] = initialConditionNextCycle(tParametricBase,nBubble,T,Y,fToSolve,Tsave,Tend,SaveHowMany,beta,Hcc,endCycle,PARAM)

    %break if there is no new nucleation happens
    breakNow = 0;
    if T(end)==Tsave(end) && nBubble==1
        disp('Cycle ended without nucleating a second bubble')
        breakNow = 1;
    elseif T(end)==Tsave(end) && nBubble==2
        disp('Cycle ended but the second bubble did not shrinks')
        breakNow = 1;
    end
    
    %compute initial condition for next cycle
    if nBubble==1   %previous simulation had one bubble, next one will have two
        
        %parameters for next iteration
        PARAM.n = [PARAM.n(1:4) PARAM.n(5) PARAM.n(5)];
        PARAM.panels = [4 1 1];
        PARAM.blockType = [1 1 1];
        PARAM.remeshProximity = {[] [] [5 6] [5 6] [3 4 6] [3 4 5]};
        
        %find postion and radius of second bubble
        Ynow = Y;
        Ynow = Ynow(end,:);
        Tnow = T;
        Tstart = Tnow(end);
        Tsave = linspace(Tstart,Tend+Tstart,SaveHowMany+1);
        posMotor = Ynow(1);
        zOld = Ynow(2);
        rOld = Ynow(3);
        [z0,r0] = findLocationAndRadiusOfSecondBubble(1,Ynow,fToSolve,beta,Hcc,PARAM);
        
        %initial conditions
        initial = [posMotor zOld z0 rOld r0]';
        
        tParametricBase(1:5) = tParametricBase(1:5);
        tParametricBase{6} = linspace(0,pi,PARAM.n(6)+1);
        
    elseif nBubble==2   %previosu simulation had two bubbles, now one or two, depending which is the end-cycle rule
        
        %from final step of previosu simulation
        Ynow = Y;
        Ynow = Ynow(end,:);
        Tnow = T;
        Tstart = Tnow(end);
        Tsave = linspace(Tstart,Tend+Tstart,SaveHowMany+1);
        
        %find out which bubble will disapper
        [~,~,~,Qmass5,Qmass6,~,~,~] = fToSolve(1,Ynow);
        
        %check if third bubble nucleates
        [nucleate,z0,r0] = findLocationAndRadiusOfThirdBubble(1,Ynow,fToSolve,beta,Hcc,PARAM);
        
        %parameters for next iteration
        if nucleate==0
            PARAM.panels = [4 1];
            PARAM.blockType = [1 1];
            PARAM.remeshProximity = {[] [] 5 5 3};
        end
        
        %if larger bubble shrinks
        if Qmass5<0 && endCycle==2
            disp('Larger bubble shrinks, eliminate it')
            
            %initial conditions
            initial = [Ynow(1) Ynow(3) Ynow(5)]';
            PARAM.n = [PARAM.n(1:4) PARAM.n(6)];
        elseif Qmass6<0 && Ynow(5)<0.3
            disp('Smaller bubble shrinks and vanishes')
            
            %initial conditions
            initial = [Ynow(1) Ynow(2) Ynow(4)]';
            PARAM.n = [PARAM.n(1:4) PARAM.n(5)];
        elseif nucleate==1
            disp('Third bubble nucleates, eliminate largest bubble')
            
            %initial conditions
            initial = [Ynow(1) Ynow(3) z0 Ynow(5) r0]';
            PARAM.n = [PARAM.n(1:4) PARAM.n(5) PARAM.n(5)];
        else
            error('Not implemented')
        end
        
        tParametricBase(1:4) = tParametricBase(1:4);
        tParametricBase{5} = linspace(0,pi,PARAM.n(5)+1);
        if nucleate==1
            tParametricBase{6} = linspace(0,pi,PARAM.n(6)+1);
        end
        
    end
    
    