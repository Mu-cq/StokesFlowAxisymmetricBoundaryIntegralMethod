%activate spectral remesh

function value = eventDragOut(t,var,flag,fVelocity)

value = 0;
if isempty(var)==0 && numel(t)==1
    
    %compute residuals
    [~,~,Un] = fVelocity(0,var);
    res = max(abs(Un));
    

    %display(['t=' num2str(t)])
    fprintf('T=%f res=%f \n',t,res);
    
else
    value = 0;
end









