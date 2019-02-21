%Runge-Kutta order 2

function [T,Y,V] = RK4(f,TTT,initial,dt,initialDT,event,PARAM)

%take variables
ttt = TTT(1):dt:TTT(end);
loop = numel(ttt);

%initialize
%loopSave = loop/PARAM.checkpoint;
Y = zeros(numel(TTT),numel(initial));
V = zeros(numel(TTT),numel(initial));
T = zeros(numel(TTT),1);

%initial condition
x1 = initial;

count = 1;
for i = 1:loop
    
    %current time
    t =ttt(i);
    
    %remesh
    if checkMetricCurvilinear(t,x1,PARAM)
        x1 = remeshDropSpectral(x1,PARAM);
    end
    
    %current initial condition
    x0 = x1;
    
    %RK4
    u0 = f(t,x0);
    u1 = f(t+0.5*dt,x0+0.5*dt*u0);
    u2 = f(t+0.5*dt,x0+0.5*dt*u1);
    u3 = f(t+dt,x0+u2*dt);
    
    %advance x0 with u1
    x1 = x0 + dt*u3;
    
    %every few iterations
    if sum(t==TTT)
        
        %staus function
        status = event(t,x0,'');
        if status==1
            break;
        end
        
        %store results
        T(count) = t;
        Y(count,:) = x0;
        V(count,:) = u0;
        
        %increase counter
        count = count+1;
        
    end

end