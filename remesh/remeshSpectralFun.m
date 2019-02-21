%remesh spectral

function [x,y,halt] = remeshSpectralFun(t,x,y,SPECTRAL)

n = numel(t);
halt = 0;

if SPECTRAL.legendre==1
    
    error('Not implemented')

%     if SPECTRAL.remeshSolve==1
%         
%         %newton method problem
%         l0 = computeTotalArcLengthSpectral(x(t),y(t),SPECTRAL);
%         fToSolve = @(tUnk) remeshSolve(tUnk,t,x,y,l0,SPECTRAL);
%         %options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','TolFun',1e-6,'TolX',1e-6,'Display','none');
%         options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','none','MaxIter',SPECTRAL.MaxIter);
%         [tNew,~,~,output] = fsolve(fToSolve,t,options);
%         if strcmp(output.message(1:17),'No solution found') && SPECTRAL.checkRemeshFail==1
%             halt = 1;
%         end
%         
%         elsec
%     
%         %constrained optimization problem
%         l0 = computeTotalArcLengthSpectral(x(t),y(t),SPECTRAL);
%         fToMinimize = @(tUnk) remeshMinimize(tUnk,t,x,y,SPECTRAL);
%         fConstr = @(tUnk) constraintArcLenght(tUnk,x,y,l0,SPECTRAL);
%         options = optimoptions('fmincon','DerivativeCheck','off','GradConstr','off','Display','notify','Algorithm','sqp','TolCon',1e-6,'TolX',1e-10,'TolFun',1e-6,'MaxIter',SPECTRAL.MaxIter,'MaxFunEvals',1e4);
%         A = []; b = []; Aeq = []; beq = []; lb = zeros(1,n); ub = ones(1,n);
%         tNew = fmincon(fToMinimize,t,A,b,Aeq,beq,lb,ub,fConstr,options);
%     
%     end  
    
elseif SPECTRAL.legendre==0||SPECTRAL.legendre==2
    
    if SPECTRAL.remeshSolve==1
        
        %newton method problem
        l0 = computeTotalArcLengthSpectral(x(t),y(t),SPECTRAL);
        fToSolve = @(tUnk) remeshSolve(tUnk,t,x,y,l0,SPECTRAL);
        options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','none','MaxIter',SPECTRAL.MaxIter);
        [tNew,~,~,output] = fsolve(fToSolve,t,options);
        if strcmp(output.message(1:17),'No solution found') && SPECTRAL.checkRemeshFail==1
            halt = 1;
        end
        %halt = 1;
        
    else
        
        %constrained optimization problem
        l0 = computeTotalArcLengthSpectral(x(t),y(t),SPECTRAL);
        fToMinimize = @(tUnk) remeshMinimize(tUnk,t,x,y,SPECTRAL);
        fConstr = @(tUnk) constraintArcLenght(tUnk,x,y,l0,SPECTRAL);
        options = optimoptions('fmincon','Display','notify','Algorithm','sqp','TolCon',1e-6,'TolX',1e-10,'TolFun',1e-6,'MaxIter',SPECTRAL.MaxIter,'MaxFunEvals',1e4);
        A = []; b = []; Aeq = []; beq = []; lb = zeros(1,n); ub = ones(1,n);
        tNew = fmincon(fToMinimize,t,A,b,Aeq,beq,lb,ub,fConstr,options);
    
    end
    
end

%new coordinates
x = x(tNew);
y = y(tNew);

%force symmetry condition
if SPECTRAL.legendre==0||SPECTRAL.legendre==2
    y([1 end]) = [0 0];
end
