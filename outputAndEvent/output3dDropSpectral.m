%activate spectral remesh

function value = output3dDropSpectral(t,var,flag,fCompute)

value = 0;
if isempty(var)==0 && numel(t)==1
    
    %check residuals
    res = 1;
    if t(1)>0 && value==0
        [~,res] = fCompute(t,var);
        res = norm(res,Inf);
    end
    
    if value==0

        fprintf('T=%f, res=%1.10f \n',t,res);

    end
    
    if res>2
        disp('Simulation crashes')
        value = 1;
    end

end









