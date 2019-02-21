%Runge-Kutta order 2

function [T,Y] = odeMatlabRemeshBubble(f,TTT,initial,event,remeshInLoop,options,ODE,PARAM)

if (TTT(2)-TTT(1))/10<options.MaxStep
   warning('Save frequency is too high compared to time stepping') 
end

%initialize
TsaveFirst = TTT;
Y = cell(numel(TTT),1);
T = zeros(numel(TTT),1);

%initial condition
x0 = initial;

count = 1;
Tvect = 0;
while Tvect(end)<TTT(end)
    
    %remesh with splines (deformable interfaces)
    [x0,remeshError] = remeshInLoop(Tvect(end),x0);

    %integrate in time
    if ODE==1
        [Tvect,Yvect] = ode45(f,TTT,x0,options);
    elseif ODE==2
        [Tvect,Yvect] = ode23t(f,TTT,x0,options);
    else
        error('Not implemented')
    end
    
    %prepare new TTT and initial condition
    x0 = Yvect(end,:)';
    [~,indMin] = min(abs(Tvect(end-1)-TTT));
    TTT = [Tvect(end) TTT(indMin+1:end)];
    
    %status function
    status = event(Tvect(end),x0,'flag');
    if status==1
            %store results
            countSave = 0;
            T(count:numel(Tvect)+count-1) = Tvect;
            for k = 1:size(Yvect,1)
                Y{count+countSave} = Yvect(k,:)';
                countSave = countSave+1;
            end
            T = T(1:numel(Tvect)+count-1);
            Y = Y(1:numel(Tvect)+count-1);
            break;
    end
        
    %store results
    countSave = 0;
    T(count:numel(Tvect)+count-1) = Tvect;
    for k = 1:size(Yvect,1)
        Y{count+countSave} = Yvect(k,:)';
        countSave = countSave+1;
    end
    
    %save results
    if PARAM.SaveDataIte==1
        
       %save results
       display('Save results')
       cd(PARAM.res)
       save(PARAM.filename)
       cd(PARAM.here)
        
    end
        
    %increase counter
    count = count+numel(Tvect);
    
    %exit when remesh add to mant points
    if remeshError==1
        break;
    end

end

%clean output data
count = 0;
for i = 1:numel(T)
   if min(abs(T(i)-TsaveFirst))>1e-5
       i = i+1;
   else
       count = count+1;
   end
   Tnew(count) = T(i);
   Ynew{count} = Y{i};
end
T = Tnew';   Y = Ynew;







