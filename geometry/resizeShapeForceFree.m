%resize shape such that irt respect force free condition

function xNew = resizeShapeForceFree(x,y,gamma)
    
%compute resizing factor
fStressX = @(unk) rescaleFunctionShape(x,y,gamma,unk,0);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','iter');
alpha = fsolve(fStressX,1,options);

%options = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolCon',1e-16,'TolX',1e-16,'TolFun',1e-16,'MaxFunEvals',1e5);

%solve non linear optimization problem subject to nonlinear constraint
%A = []; b = []; Aeq = []; beq = []; lb = 0.8; ub = 1.2;

%alpha = fmincon(fStressX,1,A,b,Aeq,beq,lb,ub);

%compute new dropet corrdinates
xNew = x*alpha;