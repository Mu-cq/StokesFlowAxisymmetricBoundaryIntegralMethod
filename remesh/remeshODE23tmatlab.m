%remesh when Using ODE from Matlab

function [T,Y] = remeshODE23tmatlab(f,Tsave,initial,options,PARAM)

TTT = Tsave;
Ysave = zeros(numel(Tsave),2*PARAM.n+2);
%loop for event
count = 0;
countRemesh = 0;
Tnow = -1;
while Tnow(end)~=Tsave(end)
    [Tnow,Ynow] = ode23t(f,TTT,initial,options);
        
    status = 0;
    %if last time is on the saving schedule
    if abs(Tnow(end)-Tsave(end))<1e-10
        display('Target time has been reached, stop simulation')
        break;
    elseif sum(Tnow(end)==Tsave)==1
        
        %save data
        if countRemesh>0
            Tnow = Tnow(2:end);
            Ynow = Ynow(2:end,:);
        end
        Ysave(count+1:count+numel(Tnow),:) = Ynow(1:end,:);
        
        %counter
        countRemesh = countRemesh + 1;
        count = count + numel(Tnow)-1;
        
        %compute residuals
        V0 = PARAM.V0;
        status = checkCoeffsAndConvergence(Tnow(end),Ynow(end,:)','',PARAM,V0);
        if status==1
            cut = find(Tsave==Tnow(end),1,'first');
            Tsave = Tsave(1:cut);
            Ysave = Ysave(1:cut,:);
        end
        
        %display('Do not execute remesh')
        
    elseif sum(Tnow(end)==Tsave)==0
        
        if countRemesh>0
            Tnow = Tnow(2:end);
            Ynow = Ynow(2:end,:);
        end
            
        ind = find(Tsave>Tnow(end),1,'first');
        TTT = [Tnow(end) Tsave(ind:end)];
            
        %save data
        Ysave(count+1:count+numel(Tnow)-1,:) = Ynow(1:end-1,:);
            
        %counter
        countRemesh = countRemesh + 1;
        count = count + numel(Tnow)-1;
       
        %execute remesh
        display(['Remesh t=' num2str(Tnow(end))])
        initial = remeshDropSpectral(Ynow(end,:)',PARAM);
        
        %compute residuals
        V0 = PARAM.V0;
        status = checkCoeffsAndConvergence(Tnow(end),Ynow(end,:)','',PARAM,V0);
        if status==1
            cut = find(Tsave>Tnow(end),1,'first')-1;
            Tsave = Tsave(1:cut);
            Ysave = Ysave(1:cut,:);
        end
        
    end
    
    if status==1
        display('Stop DNS')
        break;
    end
        
end
    
%final data
T = Tsave;
Y = Ysave;