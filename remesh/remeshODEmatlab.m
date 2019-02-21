%remesh when Using ODE from Matlab

function [T,Y] = remeshODEmatlab(f,Tsave,initial,options)

TTT = Tsave;
Ysave = zeros(numel(Tsave),2*PARAM.n+2);
%loop for event
count = 0;
Tnow = -1;
while Tnow(end)~=Tsave(end)
    [Tnow,Ynow] = ode45(f,TTT,initial,options);
        
    %if last time is on the saving schedule
    if sum(Tnow(end)==Tsave)
            
    else
            
        ind = find(Tsave>Tnow(end),1,'first');
        TTT = [Tnow(end) Tsave(ind:end)];
            
        %save data
        Ysave(count+1:count+numel(Tnow)-1,:) = Ynow(1:end-1,:);
            
        %counter
        count = count + numel(Tnow)-1;
            
    end
        
    %execute remesh
    initial = remeshDropSpectral(Ynow(end,:)',PARAM);
        
end
    
%final data
T = Tsave;
Y = Ysave;